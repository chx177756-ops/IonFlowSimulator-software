#!/bin/bash
# batch_np_scale.sh — CPU TAKUMA 扩展性测试 NP=1~10
# 每个 (mesh, NP) 跑 MAX_ITER 次迭代, 取后 AVG_START~MAX_ITER 平均 TOTAL TAKUMA Time

PDIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$PDIR/build"
RESULTS_DIR="$BUILD_DIR/batch_np_results"
MESH_DIR="$BUILD_DIR/MESH"
CONFIG_TEMPLATE="$BUILD_DIR/config.json"
MAX_ITER=10
AVG_START=6

mkdir -p "$RESULTS_DIR"
shopt -s nullglob
MESH_FILES=("$MESH_DIR"/*.msh)

if [ ${#MESH_FILES[@]} -eq 0 ]; then
    echo "ERROR: 没有网格文件 ($MESH_DIR/*.msh)"
    exit 1
fi

echo "网格数: ${#MESH_FILES[@]}, 每个跑 $MAX_ITER 次迭代, NP=1..10"
echo "取迭代 $AVG_START~$MAX_ITER 的 TOTAL TAKUMA Time 平均值"
echo "========================================"

SUMMARY="$RESULTS_DIR/summary.txt"
printf "%-35s %4s %12s\n" "网格" "NP" "Takuma(s)" > "$SUMMARY"
printf "%-35s %4s %12s\n" "---" "--" "----------" >> "$SUMMARY"

for MESH_PATH in "${MESH_FILES[@]}"; do
    MESH_NAME=$(basename "$MESH_PATH" .msh)

    sed "s|\"mesh_file\":.*|\"mesh_file\": \"./MESH/${MESH_NAME}.msh\",|" \
        "$CONFIG_TEMPLATE" > "$BUILD_DIR/config_tmp.json"
    mv "$BUILD_DIR/config_tmp.json" "$BUILD_DIR/config.json"

    for NP in $(seq 1 8); do
        echo ""
        echo ">>> [$MESH_NAME] NP=$NP"

        LOG="$RESULTS_DIR/${MESH_NAME}_np${NP}.log"
        cd "$BUILD_DIR"

        mpirun --oversubscribe -np "$NP" ./Takuma_DoubleStack_CPU 2>&1 | awk -v max="$MAX_ITER" '
            /^Iteration:/ { count++; if (count > max) exit }
            { print }
        ' > "$LOG"

        grep "TOTAL TAKUMA Time" "$LOG" > "$RESULTS_DIR/_tmp_tot.txt"
        AVG_TOT="N/A"
        if [ -s "$RESULTS_DIR/_tmp_tot.txt" ]; then
            AVG_TOT=$(tail -n +"$AVG_START" "$RESULTS_DIR/_tmp_tot.txt" | \
                      awk '{n++;sum+=$(NF-1)} END {if(n>0) printf "%.6f",sum/n; else print "N/A"}')
        fi

        printf "%-35s %4d %12s\n" "$MESH_NAME" "$NP" "$AVG_TOT" | tee -a "$SUMMARY"
        rm -f "$RESULTS_DIR"/_tmp_*.txt
    done
done

echo ""
echo "===== 结果已保存到 $SUMMARY ====="
cat "$SUMMARY"
