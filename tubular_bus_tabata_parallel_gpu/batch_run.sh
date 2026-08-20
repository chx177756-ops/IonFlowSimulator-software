#!/bin/bash
# batch_run.sh — 批量运行不同网格, 采集 Tabata 求解器性能数据
# 用法: cd /home/hxc1/projects/tubular_bus_tabata_parallel_trygpu && bash batch_run.sh

PROJECT_DIR="/home/hxc1/projects/tubular_bus_tabata_parallel_trygpu"
BUILD_DIR="$PROJECT_DIR/build"
CONFIG_TEMPLATE="$BUILD_DIR/config.json"
RESULTS_DIR="$BUILD_DIR/batch_results"
MESH_DIR="$BUILD_DIR/MESH"

MAX_ITER=10          # 每个网格跑多少迭代
AVG_START=6          # 从第几次迭代开始取平均 (迭代6~20 共15次)
NP=1                 # MPI 进程数
GPU_LIB_PATH="/home/hxc1/local/cuda-12.9/targets/x86_64-linux/lib"

mkdir -p "$RESULTS_DIR"

# ── 收集所有网格文件 ──
shopt -s nullglob
MESH_FILES=("$MESH_DIR"/*.msh)

if [ ${#MESH_FILES[@]} -eq 0 ]; then
    echo "ERROR: 没有找到网格文件 ($MESH_DIR/*.msh)"
    exit 1
fi

echo "找到 ${#MESH_FILES[@]} 个网格文件"
echo "每个跑 $MAX_ITER 次迭代, 计算 $AVG_START~$MAX_ITER 次迭代的 TOTAL Function Time 平均值"
echo "========================================"

for MESH_PATH in "${MESH_FILES[@]}"; do
    MESH_NAME=$(basename "$MESH_PATH" .msh)
    echo ""
    echo ">>> [$MESH_NAME] 开始..."

    # ── 替换 config.json 中的 mesh_file ──
    sed "s|\"mesh_file\":.*|\"mesh_file\": \"MESH/${MESH_NAME}.msh\",|" \
        "$CONFIG_TEMPLATE" > "$BUILD_DIR/config_tmp.json"
    mv "$BUILD_DIR/config_tmp.json" "$BUILD_DIR/config.json"

    # ── 运行, 第20次迭代后自动终止 ──
    LOG_FILE="$RESULTS_DIR/${MESH_NAME}.log"

    cd "$BUILD_DIR"
    export LD_LIBRARY_PATH="$GPU_LIB_PATH:$LD_LIBRARY_PATH"

    mpirun -np "$NP" ./Tubular 2>&1 | awk -v max="$MAX_ITER" '
        /^Iteration:/ { count++; if (count > max) exit }
        { print }
    ' > "$LOG_FILE"

    # ── 提取每个迭代的 TOTAL Function Time ──
    grep "TOTAL Function Time" "$LOG_FILE" > "$RESULTS_DIR/${MESH_NAME}_times.txt"

    # ── 计算第AVG_START~MAX_ITER次迭代的平均 ──
    if [ -s "$RESULTS_DIR/${MESH_NAME}_times.txt" ]; then
        AVG=$(tail -n +"$AVG_START" "$RESULTS_DIR/${MESH_NAME}_times.txt" | \
              awk '{ n++; sum += $5 } END { if(n>0) printf "%.4f", sum/n; else print "N/A" }')
        echo "$MESH_NAME: avg_total_time = $AVG s" | tee "$RESULTS_DIR/${MESH_NAME}_avg.txt"
    else
        echo "$MESH_NAME: 无 TOTAL Function Time 数据 (运行失败?)"
    fi

    echo "<<< [$MESH_NAME] 完成"
done

# ── 汇总 ──
echo ""
echo "========================================"
echo "汇总结果"
echo "========================================"
SUMMARY="$RESULTS_DIR/summary.txt"
printf "%-40s %s\n" "网格" "平均 TOTAL Time (迭代${AVG_START}~${MAX_ITER})" > "$SUMMARY"
printf "%-40s %s\n" "----" "------------------------------------" >> "$SUMMARY"
for f in "$RESULTS_DIR"/*_avg.txt; do
    name=$(basename "$f" _avg.txt)
    avg=$(cat "$f")
    printf "%-40s %s\n" "$name" "$avg" >> "$SUMMARY"
done
cat "$SUMMARY"
