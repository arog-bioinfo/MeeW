#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd -- "${script_dir}/.." && pwd -P)"
modules_root="${MEEW_MODULES_ROOT:-$(dirname -- "${repo_root}")}"

if [[ ! -d "${modules_root}" ]]; then
	printf 'Error: modules directory does not exist: %s\n' "${modules_root}" >&2
	exit 1
fi

setup_module() {
	local name="$1"
	local url="$2"
	local destination="${modules_root}/${name}"

	if [[ -e "${destination}" || -L "${destination}" ]]; then
		if [[ -d "${destination}" && -e "${destination}/.git" && -f "${destination}/workflow/Snakefile" ]]; then
			printf 'Skipping %s: valid module checkout exists at %s\n' "${name}" "${destination}"
			return
		fi

		printf 'Error: destination exists but is not a valid module checkout: %s\n' "${destination}" >&2
		printf 'Expected a directory containing .git and workflow/Snakefile; inspect it manually and remove it only if it is a partial clone.\n' >&2
		exit 1
	fi

	printf 'Cloning %s into %s\n' "${name}" "${destination}"
	git clone --branch main --single-branch -- "${url}" "${destination}"
}

setup_module "MeeW-Preprocessing" "https://github.com/arog-bioinfo/MeeW-Preprocessing.git"
setup_module "MeeW-Assembly" "https://github.com/arog-bioinfo/MeeW-Assembly.git"
setup_module "MeeW-Binning" "https://github.com/arog-bioinfo/MeeW-Binning.git"
setup_module "MeeW-Annotation" "https://github.com/arog-bioinfo/MeeW-Annotation.git"
setup_module "MeeW-Metabolic_Modeling" "https://github.com/arog-bioinfo/MeeW-Metabolic_Modeling.git"
