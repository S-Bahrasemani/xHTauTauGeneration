 #!/bin/bash

function sfu()
{
    # setup ROOT, Python etc.
    #export ROOT_VERSION=5.34.09.hsg7 #5.34.07
    # export ROOT_VERSION=5.34.10 # skims
    #export ROOT_VERSION=5.34.head
    #export ROOT_VERSION=5.34.14
    export ROOT_VERSION=5.34.18 # new TMVA
    source /atlas/software/bleedingedge/setup.sh
    export PATH=${HOME}/.local/bin${PATH:+:$PATH}

    # setup Python user base area
    export PYTHONUSERBASE=/cluster/data10/endw/local/sl5
    #export PYTHONUSERBASE=${homenodepath}/local/${OS_VERSION}
    export PATH=${PYTHONUSERBASE}/bin${PATH:+:$PATH}
}

sfu