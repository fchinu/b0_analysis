set -o pipefail
cd "$(dirname "$0")"/..
function pinfo() { echo -e "\033[32m${1}\033[m" >&2; }

test-pylint() {
    pinfo "running test: pylint"
    type pylint
    find . -name '*.py' | xargs pylint
}

test-all() {
    test-pylint
}

# Check parameters
[[ $# == 0 ]] && test-all
while [[ $# -gt 0 ]]; do
    case "$1" in

    all) test-all ;;
    pylint) test-pylint ;;
    esac
    shift
done