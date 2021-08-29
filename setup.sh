#!/bin/bash

(
set -euo pipefail

cd $(dirname ${BASH_SOURCE[0]})

PKGDIR="$(readlink -f .)"
declare -i forceStandalone=0
declare -i doPrintEnv=0
declare -i doPrintEnvInstr=0
declare -i usingCMSSW=0
declare -i needROOFITSYS_ROOTSYS=0
declare -a setupArgs=()

for farg in "$@"; do
  fargl="$(echo $farg | awk '{print tolower($0)}')"
  if [[ "$fargl" == "standalone" ]]; then
    forceStandalone=1
  elif [[ "$fargl" == "env" ]]; then
    doPrintEnv=1
  elif [[ "$fargl" == "envinstr" ]]; then
    doPrintEnvInstr=1
  else
    setupArgs+=( "$farg" ) 
  fi
done
declare -i nSetupArgs
nSetupArgs=${#setupArgs[@]}

if [[ ${forceStandalone} -eq 0 ]] && [[ ! -z "${CMSSW_BASE+x}" ]]; then

  usingCMSSW=1

  eval $(scram ru -sh)

fi

printenv() {
  # Print the environment variables from MELA as well if they are needed.
  if [[ -d ../../JHUGenMELA ]]; then
    envopts="env"
    if [[ ${forceStandalone} -eq 1 ]]; then
      envopts="${envopts} standalone"
    fi
    ../../JHUGenMELA/setup.sh ${envopts}
    eval $(../../JHUGenMELA/setup.sh ${envopts})
  fi

  if [[ ${usingCMSSW} -eq 1 ]]; then
    return 0
  fi

  libappend="${PKGDIR}CandidateLOCaster/lib:${PKGDIR}EventContainer/lib:${PKGDIR}GenericMEComputer/lib"
  end=""
  if [[ ! -z "${LD_LIBRARY_PATH+x}" ]]; then
    end=":${LD_LIBRARY_PATH}"
  fi
  if [[ "${end}" != *"$libappend"* ]]; then
    echo "export LD_LIBRARY_PATH=${libappend}${end}"
  fi
}
doenv() {
  # Set up the environment variables from MELA as well if they are needed.
  if [[ -d ../../JHUGenMELA ]]; then
    envopts="env"
    if [[ ${forceStandalone} -eq 1 ]]; then
      envopts="${envopts} standalone"
    fi
    eval $(../../JHUGenMELA/setup.sh ${envopts})
  fi

  if [[ ${usingCMSSW} -eq 1 ]]; then
    return 0
  fi
}
printenvinstr () {
  if [[ ${usingCMSSW} -eq 1 ]]; then
    return 0
  fi

  echo
  echo "remember to do"
  echo
  echo 'eval $(./setup.sh env standalone)'
  echo "or"
  echo 'eval `./setup.sh env standalone`'
  echo
  echo "if you are using a bash-related shell, or you can do"
  echo
  echo './setup.sh env standalone'
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
    if [[ ${usingCMSSW} -eq 1 ]];then
      scramv1 b "${setupArgs[@]}"
    else
      make clean
    fi

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

if [[ ${usingCMSSW} -eq 1 ]]; then
  scramv1 b "${setupArgs[@]}"

  compile_status=$?
  if [[ ${compile_status} -ne 0 ]]; then
    echo "Compilation failed with status ${compile_status}."
    exit ${compile_status}
  fi
else

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

fi

printenvinstr

)