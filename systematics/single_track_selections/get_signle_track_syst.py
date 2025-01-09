import argparse
from itertools import product
import os
from concurrent.futures import ProcessPoolExecutor
import copy
import yaml
import pandas as pd
import numpy as np
import uproot
import matplotlib.pyplot as plt

def get_all_selections(config):
    
    selections_dicts = []
    for selection_name, selection in config['selections'].items():
        selections_dicts.append([])

        # Default selection
        selections_dicts[-1].append({
            'selection_name': f"{selection_name}_default",
            'selection_string': ""
        })

        for threshold in selection['thresholds']:
            selection_string = ""
            for var_name in selection['var_names']:
                if selection['threshold_dir'] == 'lower':
                    selection_string += f"{var_name} >= {threshold} and "
                elif selection['threshold_dir'] == 'upper':
                    selection_string += f"{var_name} <= {threshold} and "
                else:
                    raise ValueError(f"Invalid threshold_dir: {selection['threshold_dir']}")
            selection_string = selection_string[:-5]
            selections_dicts[-1].append({
                'selection_name': f"{selection_name}_{threshold}",
                'selection_string': selection_string
            })

    query_dicts = []
    for selections in product(*selections_dicts):
        selection_name = ""
        query = ""
        for selection in selections:
            selection_name += f"{selection['selection_name']}_"
            if selection['selection_string'] != "":
                query += f"{selection['selection_string']} and "
        
        query_dicts.append({
            "selection_name": selection_name[:-1],
            "query": query[:-5]
        })

    return query_dicts

def create_datasets(config, query_dicts):
    df_data = pd.concat([pd.read_parquet(f) for f in config['inputs']['data']])
    df_mc = pd.concat([pd.read_parquet(f) for f in config['inputs']['mc']])

    for query_dict in query_dicts:
        output_dir = os.path.join(config['output_dir'], "data", query_dict['selection_name'])
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        if query_dict['query'] != "":
            df_data_sel = df_data.query(query_dict['query'])
            df_mc_sel = df_mc.query(query_dict['query'])
        else:
            df_data_sel = df_data
            df_mc_sel = df_mc
        df_data_sel.to_parquet(os.path.join(output_dir, 'data.parquet'))
        df_mc_sel.to_parquet(os.path.join(output_dir, 'mc.parquet'))
        del df_data_sel, df_mc_sel
    del df_data, df_mc

def create_fit_configs(config, query_dicts, fit_cfg):
    for query_dict in query_dicts:
        config_mod = copy.deepcopy(fit_cfg)
        config_mod['inputs']['data'] = [os.path.join(
            config['output_dir'],
            'data',
            query_dict['selection_name'],
            'data.parquet'
        )]
        config_mod['inputs']['mc'] = [os.path.join(
            config['output_dir'],
            'data',
            query_dict['selection_name'],
            'mc.parquet'
        )]
        config_mod['outputs']['directory'] = os.path.join(
            config['output_dir'],
            'fits',
            query_dict['selection_name']
        )
        config_mod['outputs']['suffix'] = config['raw_yields_suffix']
        config_mod['cutset_file_name'] = config['configs']['cutset']
        output_path = os.path.join(
            config['output_dir'],
            'fits',
            query_dict['selection_name'],
            f"config_fit.yml"
        )
        if not os.path.exists(os.path.join(config['output_dir'], 'fits', query_dict['selection_name'])):
            os.makedirs(os.path.join(config['output_dir'], 'fits', query_dict['selection_name']))
        with open(output_path, 'w') as f:
            yaml.dump(config_mod, f, default_flow_style=False)

def create_efficiency_configs(config, query_dicts, efficiency_cfg):
    for query_dict in query_dicts:
        config_mod = copy.deepcopy(efficiency_cfg)
        config_mod['cutset_file_name'] = config['configs']['cutset']
        config_mod['eff_frac'] = 1.
        config_mod['reco_file_names'] = [os.path.join(
            config['output_dir'],
            'data',
            query_dict['selection_name'],
            'mc.parquet'
        )]
        config_mod['output_file_name'] = os.path.join(
            config['output_dir'],
            'efficiencies',
            query_dict['selection_name'],
            config['efficiencies_file_name']
        )
        output_path = os.path.join(
            config['output_dir'],
            'efficiencies',
            query_dict['selection_name'],
            f"config_efficiency.yml"
        )
        if not os.path.exists(os.path.join(config['output_dir'], 'efficiencies', query_dict['selection_name'])):
            os.makedirs(os.path.join(config['output_dir'], 'efficiencies', query_dict['selection_name']))

        with open(output_path, 'w') as f:
            yaml.dump(config_mod, f, default_flow_style=False)

def create_cross_section_configs(config, query_dicts, cross_section_cfg):
    for query_dict in query_dicts:
        config_mod = copy.deepcopy(cross_section_cfg)
        rawyields_file_name = f"{config['particle']}_mass{config['raw_yields_suffix']}.root"
        config_mod['rawyield_file'] = os.path.join(
            config['output_dir'],
            'fits',
            query_dict['selection_name'],
            rawyields_file_name
        )
        config_mod['efficiency_file'] = os.path.join(
            config['output_dir'],
            'efficiencies',
            query_dict['selection_name'],
            config['efficiencies_file_name']
        )
        config_mod['output_file'] = os.path.join(
            config['output_dir'],
            'cross_sections',
            query_dict['selection_name'],
            config["cross_section_file_name"]
        )
        output_path = os.path.join(
            config['output_dir'],
            'cross_sections',
            query_dict['selection_name'],
            f"config_cross_section.yml"
        )
        if not os.path.exists(os.path.join(config['output_dir'], 'cross_sections', query_dict['selection_name'])):
            os.makedirs(os.path.join(config['output_dir'], 'cross_sections', query_dict['selection_name']))
    
        with open(output_path, 'w') as f:
            yaml.dump(config_mod, f, default_flow_style=False)

def extract_rawyield(fit_config_name):
    os.system(f"python3 fit/extract_rawyield.py -c {fit_config_name}")

def extract_efficiency(efficiency_config_name):
    os.system(f"python3 efficiency/get_efficiency_bmesons.py {efficiency_config_name}")

def extract_cross_section(cross_section_config_name):
    os.system(f"python3 cross_section/compute_cross_section.py {cross_section_config_name}")

def draw_results(config, query_dicts):
    rawyields, efficiencies, cross_sections = [], [], []
    rawyields_unc, efficiencies_unc, cross_sections_unc = [], [], []
    fig, axs = plt.subplots(2, 2, figsize=(12, 9))
    fig.subplots_adjust(wspace = 0.02)
    fig_ratio, axs_ratio = plt.subplots(1, 2, figsize=(12, 5))
    fig_ratio.subplots_adjust(wspace = 0.02)
    colors = plt.cm.tab20.colors
    for i, query_dict in enumerate(query_dicts):
        rawyields_file_name = f"{config['particle']}_mass{config['raw_yields_suffix']}.root"
        rawy_path = os.path.join(
            config['output_dir'],
            'fits',
            query_dict['selection_name'],
            rawyields_file_name
        )
        with uproot.open(rawy_path) as f:
            if "h_rawyields" not in f:
                continue
            rawyields.append(f["h_rawyields"].values())
            rawyields_unc.append(f["h_rawyields"].errors())
            if i == 0:
                pt_bins = f["h_rawyields"].axis().edges()

        eff_path = os.path.join(    
            config['output_dir'],
            'efficiencies',
            query_dict['selection_name'],
            config['efficiencies_file_name']
        )
        with uproot.open(eff_path) as f:
            efficiencies.append(f["h_eff"].values())
            efficiencies_unc.append(f["h_eff"].errors())

        cross_section_path = os.path.join(
            config['output_dir'],
            'cross_sections',
            query_dict['selection_name'],
            config["cross_section_file_name"]
        )
        with uproot.open(cross_section_path) as f:
            cross_sections.append(f["h_cross_section"].values())
            cross_sections_unc.append(f["h_cross_section"].errors())
        
    pt_centres = (pt_bins[:-1] + pt_bins[1:])/2
    pt_widths = (pt_bins[1:] - pt_bins[:-1])/2
    ax_titles = ['Raw yields', 'Efficiency', 'Cross section', 'Counts']
    for i_trial, (rawy, eff, cross_sec, rawy_unc, eff_unc, cross_sec_unc, query_dict) in enumerate(zip(rawyields, efficiencies, cross_sections, rawyields_unc, efficiencies_unc, cross_sections_unc, query_dicts)):
        axs[0, 0].errorbar(
            pt_centres, rawy, yerr=rawy_unc, xerr=pt_widths,
            fmt='o', label=query_dict['selection_name'], c=colors[i_trial]
        )
        axs[0, 1].errorbar(
            pt_centres, eff, yerr=eff_unc, xerr=pt_widths,
            fmt='o', label=query_dict['selection_name'], c=colors[i_trial]
        )
        axs[1, 0].errorbar(
            pt_centres, cross_sec, yerr=cross_sec_unc, xerr=pt_widths,
            fmt='o', label=query_dict['selection_name'], c=colors[i_trial]
        )
        axs_ratio[0].errorbar(
            pt_centres, cross_sec, yerr=cross_sec_unc, xerr=pt_widths,
            fmt='o', label=query_dict['selection_name'], c=colors[i_trial]
        )

        ratio_unc = cross_sec/cross_sections[0] * np.sqrt((cross_sec_unc/cross_sec)**2 + (cross_sections_unc[0]/cross_sections[0])**2)
        axs_ratio[1].errorbar(
            pt_centres, cross_sec/cross_sections[0], yerr=ratio_unc, xerr=pt_widths,
            fmt='o', label=query_dict['selection_name'], c=colors[i_trial]
        )
    axs[1, 0].set_yscale('log')
    axs_ratio[0].set_yscale('log')
    axs_ratio[1].set_ylim(0.7, 1.3)
    axs_ratio[1].axhline(1, color='black', linestyle='--')

    rms_shifts = []
    for i_pt in range(len(cross_sections[0])):
        cross_sections_pt = [cross_sec[i_pt] for cross_sec in cross_sections]
        rms_shift = np.array(cross_sections_pt).std()**2
        rms_shift += (np.mean(cross_sections_pt) - cross_sections_pt[0])**2
        rms_shift = np.sqrt(rms_shift)
        rms_shifts.append(rms_shift)

    axs[1, 1].errorbar(
            pt_centres,
            cross_sections[0],
            yerr=rms_shifts,
            xerr=pt_widths,
            label=r'$\mathrm{\sqrt{RMS^2 + \Delta^2}}$',
            fmt='o'
        )
    axs[1, 1].bar(pt_centres, np.ones(pt_centres.shape) * 2 * rms_shifts, bottom=np.array(cross_sections[0])-rms_shifts, width=pt_widths, color='C0', alpha=0.3, zorder=0, align='center')
    axs[1, 1].errorbar(
            pt_centres,
            cross_sections[0],
            yerr=np.array(config['assigned_syst'])*cross_sections[0],
            xerr=pt_widths,
            label='Assigned syst.',
            fmt='o'
        )
    axs[1, 1].legend()
    axs[1, 1].set_yscale('log')

    for i_ax, ax in enumerate(axs.flat):
        ax.set_xlabel(r'$p_\mathrm{T}  (\mathrm{GeV}/c)$')
        ax.set_ylabel(ax_titles[i_ax])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    axs[0, 1].legend(loc='lower right', bbox_to_anchor=(2, -0.5), fontsize=8)
    fig.savefig(os.path.join(config['output_dir'], 'results.pdf'), bbox_inches='tight')
    
    for i_ax, ax in enumerate(axs_ratio.flat):
        ax.set_xlabel(r'$p_\mathrm{T}  (\mathrm{GeV}/c)$')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    axs_ratio[0].set_ylabel('Cross section')
    axs_ratio[1].set_ylabel('Ratio')
    
    axs_ratio[1].legend(loc='lower right', bbox_to_anchor=(2, 0.25), fontsize=8)
    fig_ratio.savefig(os.path.join(config['output_dir'], 'results_ratio.pdf'), bbox_inches='tight')

    dfs = [[] for _ in range(len(rawyields[0]))]
    for rawy, eff, cross_sec, rawy_unc, eff_unc, cross_sec_unc, query_dict in zip(rawyields, efficiencies, cross_sections, rawyields_unc, efficiencies_unc, cross_sections_unc, query_dicts):
        for i_pt in range(len(rawy)):
            dfs[i_pt].append([
                rawy[i_pt],
                rawy_unc[i_pt],
                eff[i_pt],
                eff_unc[i_pt],
                cross_sec[i_pt],
                cross_sec_unc[i_pt],
                query_dict['query']
            ])
    cols_names = ["raw_yield", "raw_yield_unc", "efficiency", "efficiency_unc", "cross_section", "cross_section_unc", "query"]
    dfs = [pd.DataFrame(df, columns=cols_names) for df in dfs]
    for df, pt_min, pt_max in zip(dfs, pt_bins[:-1], pt_bins[1:]):
        df.to_parquet(os.path.join(config['output_dir'], f'results_{pt_min*10:.0f}_{pt_max*10:.0f}.parquet'))

    with uproot.recreate(os.path.join(config['output_dir'], 'results.root')) as f:
        f['rms_shifts_sum_quadrature'] = (np.array(rms_shifts)/np.array(cross_sections[0]), np.array(pt_bins))
        f['assigned_syst'] = (np.array(config['assigned_syst']), np.array(pt_bins))


def produce_files_for_syst(config_file, data, raw_yields, efficiency, cross_section, draw):
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    with open(config["configs"]["cutset"], "r") as f:
        cutset = yaml.safe_load(f)

    with open(config["configs"]["fit"], "r") as f:
        fit_cfg = yaml.safe_load(f)

    with open(config["configs"]["efficiency"], "r") as f:
        efficiency_cfg = yaml.safe_load(f)

    with open(config["configs"]["cross_section"], "r") as f:
        cross_section_cfg = yaml.safe_load(f)

    if not os.path.exists(config['output_dir']):
        os.makedirs(config['output_dir'])

    query_dicts = get_all_selections(config)

    if data:
        create_datasets(config, query_dicts)
    if raw_yields:
        create_fit_configs(config, query_dicts, fit_cfg)
    if efficiency:
        create_efficiency_configs(config, query_dicts, efficiency_cfg)
    if cross_section:
        create_cross_section_configs(config, query_dicts, cross_section_cfg)

    if raw_yields:
        with ProcessPoolExecutor(max_workers=config["max_workers"]) as executor:
            for query_dict in query_dicts:
                fit_config_name = os.path.join(
                    config['output_dir'],
                    'fits',
                    query_dict['selection_name'],
                    "config_fit.yml"
                )
                executor.submit(extract_rawyield, fit_config_name)

    if efficiency:
        with ProcessPoolExecutor(max_workers=config["max_workers"]) as executor:
            for query_dict in query_dicts:
                efficiency_config_name = os.path.join(
                    config['output_dir'],
                    'efficiencies',
                    query_dict['selection_name'],
                    "config_efficiency.yml"
                )
                executor.submit(extract_efficiency, efficiency_config_name)
            
    if cross_section:
        with ProcessPoolExecutor(max_workers=config["max_workers"]) as executor:
            for query_dict in query_dicts:
                cross_section_config_name = os.path.join(
                    config['output_dir'],
                    'cross_sections',
                    query_dict['selection_name'],
                    "config_cross_section.yml"
                )
                executor.submit(extract_cross_section, cross_section_config_name)

    if draw:
        draw_results(config, query_dicts)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Produce files for single track systematics')
    parser.add_argument("--config", "-c", metavar="text", default="config_fit.yml",
                        help="yaml config file for fit", required=True)
    parser.add_argument("--data", help="Produce datasets", action="store_true")
    parser.add_argument("--raw-yields", "-r", help="Produce raw yields", action="store_true")
    parser.add_argument("--efficiency", "-e", help="Produce efficiencies", action="store_true")
    parser.add_argument("--cross-section", "-x", help="Produce cross sections", action="store_true")
    parser.add_argument("--draw", help="Draw results", action="store_true")
    args = parser.parse_args()

    produce_files_for_syst(
        args.config, args.data, args.raw_yields,
        args.efficiency, args.cross_section, args.draw)