#!/usr/bin/env bash
set -e

if [ $# -ne 4 ]; then
    echo "Usage: $0 <CONFIG_LIST> <TIMESTEP_LIST> <LOCAL_SIM_DIR> <REMOTE_SIM_DIR>"
    exit 1
fi

read -r -a CONFIG_LIST <<< "$1"
read -r -a TIMESTEP_LIST <<< "$2"
LOCAL_SIM_DIR="$3"
REMOTE_SIM_DIR="$4"

REMOTE_USER="kersweb"
REMOTE_HOST="barklalogin1.liv.ac.uk"
CONTROL_SOCKET="$HOME/.ssh/cm_socket_%h_%p_%r"
SSH_BASE=(ssh -o ControlMaster=auto -o ControlPersist=15m -o ControlPath="$CONTROL_SOCKET")
RSYNC_OPTS=(-e "${SSH_BASE[*]}" -azP)

echo "    --------------------------------------------------"
echo "    Syncing .run files"
echo "    --------------------------------------------------"

echo "    Establishing master SSH connection (expect Duo 2FA prompt) ..."
"${SSH_BASE[@]}" -MNf "${REMOTE_USER}@${REMOTE_HOST}"

REMOTE_CONFIGS_DIR="${REMOTE_SIM_DIR}/configs"
LOCAL_CONFIGS_DIR="${LOCAL_SIM_DIR}/configs"

mkdir -p "$LOCAL_CONFIGS_DIR"
rsync "${RSYNC_OPTS[@]}" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_CONFIGS_DIR}/*.run" "$LOCAL_CONFIGS_DIR/" 2>/dev/null || echo " !! Warning: no run files found"

for prm_path in "${CONFIG_LIST[@]}"; do
    model_name="$(basename "$prm_path" .prm | sed -E 's/-/_/g; s/(lnk)_([+-]?[0-9]+)/\1-\2/g')"

    echo "    --------------------------------------------------"
    echo "    Syncing model: $model_name"
    echo "    --------------------------------------------------"

    REMOTE_RESULTS_DIR="${REMOTE_SIM_DIR}/results/$model_name"
    LOCAL_RESULTS_DIR="${LOCAL_SIM_DIR}/results/$model_name"

    mkdir -p "$LOCAL_RESULTS_DIR/solution"

    rsync "${RSYNC_OPTS[@]}" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/statistics" "$LOCAL_RESULTS_DIR/" \
      2>/dev/null || echo " !! Warning: no statistics file found"
    rsync "${RSYNC_OPTS[@]}" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/depth_average.txt" "$LOCAL_RESULTS_DIR/" \
      2>/dev/null || echo " !! Warning: no depth_average.txt file found"

    for tstep in "${TIMESTEP_LIST[@]}"; do
        rsync "${RSYNC_OPTS[@]}" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/solution/solution-00*${tstep}.*.vtu" "$LOCAL_RESULTS_DIR/solution/" \
          2>/dev/null || echo " !! Warning: no vtu solution files for timestep $tstep"
        rsync "${RSYNC_OPTS[@]}" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/solution/solution-00*${tstep}.pvtu" "$LOCAL_RESULTS_DIR/solution/" \
          2>/dev/null || echo " !! Warning: no pvtu solution files for timestep $tstep"
    done
done

echo "    --------------------------------------------------"
echo "    All syncs complete!"
