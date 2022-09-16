#!/bin/bash

(
set -euo pipefail

cd $(dirname ${BASH_SOURCE[0]})

PKGDIR="$(readlink -f .)"
declare -i doPrintEnv=0
declare -i doPrintEnvInstr=0
declare -i needROOFITSYS_ROOTSYS=0
declare -a setupArgs=()

for farg in "$@"; do
  fargl="$(echo $farg | awk '{print tolower($0)}')"
  if [[ "$fargl" == "env" ]]; then
    doPrintEnv=1
  elif [[ "$fargl" == "envinstr" ]]; then
    doPrintEnvInstr=1
  else
    setupArgs+=( "$farg" ) 
  fi
done
declare -i nSetupArgs
nSetupArgs=${#setupArgs[@]}

if [[ ! -d ../JHUGenMELA ]]; then
  echo "${PKGDIR}/../JHUGenMELA" does not exist.
  exit 1
fi

printenv() {
  ../JHUGenMELA/setup.sh env
  eval $(../JHUGenMELA/setup.sh env)

  if [[ -z "${MELAANALYTICS_PATH+x}" ]] || [[ "${MELAANALYTICS_PATH}" != "${PKGDIR}" ]]; then
    echo "export MELAANALYTICS_PATH=${PKGDIR}"
    export MELAANALYTICS_PATH="${PKGDIR}"
  fi

  libappend="${MELAANALYTICS_PATH}/CandidateLOCaster/lib:${MELAANALYTICS_PATH}/EventContainer/lib:${MELAANALYTICS_PATH}/GenericMEComputer/lib"
  end=""
  if [[ ! -z "${LD_LIBRARY_PATH+x}" ]]; then
    end=":${LD_LIBRARY_PATH}"
  fi
  if [[ "${end}" != *"$libappend"* ]]; then
    echo "export LD_LIBRARY_PATH=${libappend}${end}"
  fi
}
doenv() {
  eval $(../JHUGenMELA/setup.sh env)

  if [[ -z "${MELAANALYTICS_PATH+x}" ]] || [[ "${MELAANALYTICS_PATH}" != "${PKGDIR}" ]]; then
    export MELAANALYTICS_PATH="${PKGDIR}"
  fi
}
printenvinstr () {
  echo
  echo "remember to do"
  echo
  echo 'eval $('${BASH_SOURCE[0]}' env)'
  echo "or"
  echo 'eval `'${BASH_SOURCE[0]}' env`'
  echo
  echo "if you are using a bash-related shell, or you can do"
  echo
  echo ${BASH_SOURCE[0]}' env'
  echo
  echo "and change the commands according to your shell in order to do something equivalent to set up the environment variables."
  echo
}

if [[ $doPrintEnv -eq 1 ]]; then
    printenv
    exit
elif [[ $doPrintEnvInstr -eq 1 ]]; then
    printenvinstr
    exit
fi

if [[ $nSetupArgs -eq 0 ]]; then
    setupArgs+=( -j 1 )
    nSetupArgs=2
fi


if [[ "$nSetupArgs" -eq 1 ]] && [[ "${setupArgs[0]}" == *"clean"* ]]; then
    for ff in $(find ./ -name makefile); do
      ff=${ff//'/makefile'}
      cd $ff &> /dev/null
      make clean
      cd - &> /dev/null
    done

    exit $?
elif [[ "$nSetupArgs" -ge 1 ]] && [[ "$nSetupArgs" -le 2 ]] && [[ "${setupArgs[0]}" == *"-j"* ]]; then
    : ok
else
    echo "Unknown arguments:"
    echo "  ${setupArgs[@]}"
    echo "Should be nothing, env, or clean"
    exit 1
fi

doenv

for ff in $(find ./ -name makefile); do
  ff=${ff//'/makefile'}
  cd $ff &> /dev/null

  make "${setupArgs[@]}"

  compile_status=$?
  if [[ ${compile_status} -ne 0 ]]; then
    echo "Compilation of ${ff} failed with status ${compile_status}."
    exit ${compile_status}
  fi

  cd - &> /dev/null
done

printenvinstr

)
