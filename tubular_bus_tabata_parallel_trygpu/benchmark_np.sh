#!/bin/bash
# 并行分区数对加载耗时的研究 (同一网格, 不同进程数)
# 目录结构: 进程数量研究/<MeshSize>/<NP>/mesh<MeshSize>_*.msh
# 用法: bash tools/benchmark_np.sh

set -e

BASE_DIR="/home/hxc1/projects/tubular_bus_tabata_parallel_trygpu/build"
MESH_DIR="$BASE_DIR/MESH/进程数量研究"
CONVERTER="$BASE_DIR/gmsh_partition_to_mfem"
TUBULAR="$BASE_DIR/Tubular"
RESULT_FILE="$BASE_DIR/benchmark_np_results.txt"
REPEATS=5

echo "=============================================="
echo "  并行分区数对加载耗时影响"
echo "  重复次数: $REPEATS"
echo "=============================================="
echo ""

printf "%-10s %-3s %-4s %-10s %-10s\n" "MeshSize" "NP" "Run" "ConvTime_s" "LoadTime_s" | tee "$RESULT_FILE"
echo "----------------------------------------------------------------" | tee -a "$RESULT_FILE"

for mesh_dir in "$MESH_DIR"/*/; do
    mesh=$(basename "$mesh_dir")
    [[ "$mesh" =~ ^[0-9]+$ ]] || continue

    for np_dir in "$mesh_dir"*/; do
        np=$(basename "$np_dir")
        [[ "$np" =~ ^[0-9]+$ ]] || continue

        # 收集 GMSH 文件
        if [ "$np" -eq 1 ]; then
            files=$(ls "$np_dir"*.msh 2>/dev/null | grep -v "^\._")
            nfiles=$(echo "$files" | wc -w)
        else
            files=$(ls "$np_dir"mesh${mesh}_*.msh 2>/dev/null | grep -v "^\._" | sort -V)
            nfiles=$(echo "$files" | wc -w)
        fi

        if [ "$nfiles" -ne "$np" ]; then
            echo "  [跳过] Mesh=$mesh NP=$np: 文件数=$nfiles, 期望=$np"
            continue
        fi

        prefix="$np_dir/mesh${mesh}.msh.part"

        echo ""
        echo "=== Mesh=$mesh, NP=$np ==="

        for run in $(seq 1 $REPEATS); do
            # ---- 1. 转换耗时 (np=1 时跳过) ----
            if [ "$np" -eq 1 ]; then
                t_conv="0.0000"
            else
                output=$(mpirun --oversubscribe -np $np $CONVERTER $files -o $prefix 2>&1)
                t_conv=$(echo "$output" | grep "Total time" | awk '{print $3}')
            fi

            # ---- 2. 加载耗时 ----
            CONFIG_TMP="/tmp/config_bench_$$.json"
            cp "$BASE_DIR/config.json" "$CONFIG_TMP"
            python3 -c "
import json
with open('$CONFIG_TMP') as f:
    c = json.load(f)
c['mesh_file'] = '${np_dir}/mesh${mesh}.msh'
with open('$CONFIG_TMP', 'w') as f:
    json.dump(c, f, indent=4)
" 2>/dev/null

            output2=$(MESH_BENCH=1 mpirun --oversubscribe -np $np $TUBULAR -c "$CONFIG_TMP" 2>&1)
            t_load=$(echo "$output2" | grep "MESH_LOAD_TIME" | awk '{print $2}')
            rm -f "$CONFIG_TMP"

            printf "  Run %d: conv=%ss load=%ss\n" $run "$t_conv" "$t_load"
            printf "%-10d %-3d %-4d %-10s %-10s\n" $mesh $np $run "$t_conv" "$t_load" | tee -a "$RESULT_FILE"
        done
    done
done

echo ""
echo "=============================================="
echo "  测试完成。结果: $RESULT_FILE"
echo "=============================================="

# 平均值汇总
python3 - "$RESULT_FILE" << 'PYEOF'
import sys
results = {}
with open(sys.argv[1]) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) != 5: continue
        mesh, np, run, tconv, tload = parts
        try:
            key = (int(mesh), int(np))
            results.setdefault(key, {'conv': [], 'load': []})
            if tconv != '0.0000':
                results[key]['conv'].append(float(tconv))
            results[key]['load'].append(float(tload))
        except: pass

print("\n========== 平均值汇总 ==========")
print(f"{'Mesh':<10} {'NP':<4} {'Conv_avg(s)':<12} {'Load_avg(s)':<12}")
for (mesh, np), data in sorted(results.items()):
    ca = sum(data['conv'])/len(data['conv']) if data['conv'] else 0
    la = sum(data['load'])/len(data['load'])
    conv_str = f"{ca:<12.4f}" if data['conv'] else "     N/A    "
    print(f"{mesh:<10} {np:<4} {conv_str} {la:<12.4f}")
PYEOF "$RESULT_FILE"
