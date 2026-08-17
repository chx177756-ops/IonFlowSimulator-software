#!/bin/bash
# batch_np_scale.sh — 测试 NP=1~10 的 GPU Tabata 扩展性
# 每个 (mesh, NP) 跑 10 次迭代, 取后 5 次平均值
# 指标: Assembly Loop, GMRES+AMG(Setup+Solve), TOTAL Function Time

PDIR="/home/hxc1/projects/tubular_bus_tabata_parallel_trygpu"
BUILD_DIR="$PDIR/build"
RESULTS_DIR="$BUILD_DIR/batch_np_results"
MESH_DIR="$BUILD_DIR/MESH"
CONFIG_TEMPLATE="$BUILD_DIR/config.json"
MAX_ITER=5
AVG_START=2
GPU_LIB_PATH="/home/hxc1/local/cuda-12.9/targets/x86_64-linux/lib"

mkdir -p "$RESULTS_DIR"
shopt -s nullglob
MESH_FILES=("$MESH_DIR"/*.msh)

if [ ${#MESH_FILES[@]} -eq 0 ]; then
    echo "ERROR: 没有网格文件 ($MESH_DIR/*.msh)"
    exit 1
fi

echo "网格数: ${#MESH_FILES[@]}, 每个跑 $MAX_ITER 次迭代, NP=1..10"
echo "取迭代 $AVG_START~$MAX_ITER 的平均值"
echo "========================================"

SUMMARY="$RESULTS_DIR/summary.txt"
printf "%-35s %4s %12s %12s %12s\n" "网格" "NP" "Assembly(s)" "GMRES+AMG(s)" "Total(s)" > "$SUMMARY"
printf "%-35s %4s %12s %12s %12s\n" "---" "--" "----------" "------------" "-------" >> "$SUMMARY"

for MESH_PATH in "${MESH_FILES[@]}"; do
    MESH_NAME=$(basename "$MESH_PATH" .msh)

    sed "s|\"mesh_file\":.*|\"mesh_file\": \"MESH/${MESH_NAME}.msh\",|" \
        "$CONFIG_TEMPLATE" > "$BUILD_DIR/config_tmp.json"
    mv "$BUILD_DIR/config_tmp.json" "$BUILD_DIR/config.json"

    for NP in $(seq 1 5); do
        echo ""
        echo ">>> [$MESH_NAME] NP=$NP"

        LOG="$RESULTS_DIR/${MESH_NAME}_np${NP}.log"
        cd "$BUILD_DIR"
        export LD_LIBRARY_PATH="$GPU_LIB_PATH:$LD_LIBRARY_PATH"

        mpirun --oversubscribe -np "$NP" ./Tubular 2>&1 | awk -v max="$MAX_ITER" '
            /^Iteration:/ { count++; if (count > max) exit }
            { print }
        ' > "$LOG"

        # Assembly Loop
        grep "Matrix Assembly Loop" "$LOG" > "$RESULTS_DIR/_tmp_asm.txt"
        AVG_ASM="N/A"
        if [ -s "$RESULTS_DIR/_tmp_asm.txt" ]; then
            AVG_ASM=$(tail -n +"$AVG_START" "$RESULTS_DIR/_tmp_asm.txt" | \
                      awk '{n++; sum+=$(NF-1)} END {if(n>0) printf "%.6f",sum/n; else print "N/A"}')
        fi

        # GMRES+AMG = Setup + Solve
        grep "GMRES+AMG Setup\|GMRES+AMG Solve" "$LOG" > "$RESULTS_DIR/_tmp_gm.txt"
        AVG_GM="N/A"
        if [ -s "$RESULTS_DIR/_tmp_gm.txt" ]; then
            TL=$((AVG_START * 2 - 1))
            AVG_GM=$(tail -n +"$TL" "$RESULTS_DIR/_tmp_gm.txt" | \
                     awk '{v=$(NF-1); if(NR%2==1)s=v; else{n++;sum+=(s+v)}} END {if(n>0)printf "%.6f",sum/n; else print "N/A"}')
        fi

        # TOTAL Function Time
        grep "TOTAL Function Time" "$LOG" > "$RESULTS_DIR/_tmp_tot.txt"
        AVG_TOT="N/A"
        if [ -s "$RESULTS_DIR/_tmp_tot.txt" ]; then
            AVG_TOT=$(tail -n +"$AVG_START" "$RESULTS_DIR/_tmp_tot.txt" | \
                      awk '{n++;sum+=$(NF-1)} END {if(n>0) printf "%.6f",sum/n; else print "N/A"}')
        fi

        printf "%-35s %4d %12s %12s %12s\n" "$MESH_NAME" "$NP" "$AVG_ASM" "$AVG_GM" "$AVG_TOT" | tee -a "$SUMMARY"
        rm -f "$RESULTS_DIR"/_tmp_*.txt
    done
done

echo ""
echo "===== 结果已保存到 $SUMMARY ====="
cat "$SUMMARY"
