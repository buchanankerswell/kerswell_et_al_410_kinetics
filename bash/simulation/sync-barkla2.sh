#!/usr/bin/env bash
set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 <CONFIG_LIST> <TIMESTEP_LIST>"
    exit 1
fi

read -r -a CONFIG_LIST <<< "$1"
read -r -a TIMESTEP_LIST <<< "$2"

REMOTE_USER="kersweb"
REMOTE_HOST="barklalogin1.liv.ac.uk"
CONTROL_SOCKET="$HOME/.ssh/cm_socket_%h_%p_%r"
SSH_BASE=(ssh -o ControlMaster=auto -o ControlPersist=15m -o ControlPath="$CONTROL_SOCKET")
RSYNC_OPTS=(-e "${SSH_BASE[*]}" -azP)

echo "Establishing master SSH connection (expect Duo 2FA prompt) ..."
"${SSH_BASE[@]}" -MNf "${REMOTE_USER}@${REMOTE_HOST}"

for prm_path in "${CONFIG_LIST[@]}"; do
    model_name="$(basename "$prm_path" .prm | tr '-' '_')"

    REMOTE_CONFIGS_DIR="/users/kersweb/scratch/kerswell_et_al_dynp/simulation/2d_box/configs"
    REMOTE_RESULTS_DIR="/users/kersweb/scratch/kerswell_et_al_dynp/simulation/2d_box/results/$model_name"
    LOCAL_CONFIGS_DIR="$HOME/Working/kerswell_et_al_dynp/simulation/2d_box/configs"
    LOCAL_RESULTS_DIR="$HOME/Working/kerswell_et_al_dynp/simulation/2d_box/results/$model_name"

    echo "---------------------------------------------------------"
    echo "Syncing for model: $model_name"
    echo "---------------------------------------------------------"

    mkdir -p "$LOCAL_CONFIGS_DIR" "$LOCAL_RESULTS_DIR/solution"

    echo "Trying to rsync .out files from: ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/*.out"
    rsync "${RSYNC_OPTS[@]}" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/*.out" \
      "$LOCAL_RESULTS_DIR/" 2>/dev/null || echo "Warning: No .out files found"

    echo "Trying to rsync .run files from: ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_CONFIGS_DIR}/*.run"
    rsync "${RSYNC_OPTS[@]}" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_CONFIGS_DIR}/*.run" \
      "$LOCAL_CONFIGS_DIR/" 2>/dev/null || echo "Warning: No .run files found"

    echo "Trying to rsync 'statistics' from: ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/statistics"
    rsync "${RSYNC_OPTS[@]}" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/statistics" \
      "$LOCAL_RESULTS_DIR/" 2>/dev/null || echo "Warning: No 'statistics' file found"

    echo "Trying to rsync 'depth_average.txt' from: ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/depth_average.txt"
    rsync "${RSYNC_OPTS[@]}" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/depth_average.txt" \
      "$LOCAL_RESULTS_DIR/" 2>/dev/null || echo "Warning: 'No depth_average.txt' file found"

    echo "Trying to rsync 'stokes_residuals.txt' from: ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/stokes_residuals.txt"
    rsync "${RSYNC_OPTS[@]}" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/stokes_residuals.txt" \
      "$LOCAL_RESULTS_DIR/" 2>/dev/null || echo "Warning: No 'stokes_residuals.txt' file found"

    for tstep in "${TIMESTEP_LIST[@]}"; do
        echo "Trying to rsync solution files for timestep $tstep from: ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/solution/solution-*${tstep}*"
        rsync "${RSYNC_OPTS[@]}" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/solution/solution-*${tstep}*" \
          "$LOCAL_RESULTS_DIR/solution/" 2>/dev/null || echo "Warning: No solution files for timestep $tstep"
    done
done

echo "========================================================="
echo "All syncs complete!"
