#!/bin/bash
# batch_run.sh — CPU TAKUMA 批量运行 (awk 杀进程, 不修改 Main)
# 前置: CVPDE_TAKUMA.cpp 末尾有 TOTAL TAKUMA Time 计时输出

PDIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$PDIR/build"
RESULTS_DIR="$BUILD_DIR/batch_results"
MESH_DIR="$BUILD_DIR/MESH"
CONFIG_TEMPLATE="$BUILD_DIR/config.json"
MAX_ITER=10
AVG_START=6
NP=1

mkdir -p "$RESULTS_DIR"
shopt -s nullglob
MESH_FILES=("$MESH_DIR"/*.msh)

echo "网格数: ${#MESH_FILES[@]}, 各跑 $MAX_ITER 次迭代, NP=$NP"
echo "取迭代 $AVG_START~$MAX_ITER 的 TOTAL TAKUMA Time 平均值"
echo "========================================"

for MESH_PATH in "${MESH_FILES[@]}"; do
    MESH_NAME=$(basename "$MESH_PATH" .msh)
    echo ""
    echo ">>> [$MESH_NAME]"

    sed "s|\"mesh_file\":.*|\"mesh_file\": \"./MESH/${MESH_NAME}.msh\",|" \
        "$CONFIG_TEMPLATE" > "$BUILD_DIR/config_tmp.json"
    mv "$BUILD_DIR/config_tmp.json" "$BUILD_DIR/config.json"

    LOG="$RESULTS_DIR/${MESH_NAME}.log"
    cd "$BUILD_DIR"

    mpirun -np "$NP" ./Takuma_DoubleStack_CPU 2>&1 | awk -v max="$MAX_ITER" '
        /^Iteration:/ { count++; if (count > max) exit }
        { print }
    ' > "$LOG"

    grep "TOTAL TAKUMA Time" "$LOG" > "$RESULTS_DIR/${MESH_NAME}_times.txt"

    if [ -s "$RESULTS_DIR/${MESH_NAME}_times.txt" ]; then
        AVG=$(tail -n +"$AVG_START" "$RESULTS_DIR/${MESH_NAME}_times.txt" | \
              awk '{ n++; sum += $5 } END { if(n>0) printf "%.6f", sum/n; else print "N/A" }')
        echo "$MESH_NAME  avg_time=$AVG s" | tee "$RESULTS_DIR/${MESH_NAME}_avg.txt"
    else
        echo "$MESH_NAME  数据缺失" | tee "$RESULTS_DIR/${MESH_NAME}_avg.txt"
    fi
    echo "<<< [$MESH_NAME]"
done

echo ""
echo "===== 汇总 ====="
for f in "$RESULTS_DIR"/*_avg.txt; do cat "$f"; done | tee "$RESULTS_DIR/summary.txt"
