#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 13:23:42 2026

@author: han
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 13:11:58 2026

@author: han
"""

import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import json
import os
import matplotlib.ticker as ticker
import matplotlib.patches as patches

# 把 CSV 中的“起点—终点—作用值”关系转换成 GLMY 输入格式 + Python 带权有向图格式
# def generation_unfiltered_data(csv_path, weight_shift=100):
#     df = pd.read_csv(csv_path)
#
#     vertices = sorted(set(df['From']).union(set(df['To'])))
#     V = {name: idx + 1 for idx, name in enumerate(vertices)}
#
#     input_V = ','.join(str(V[name]) for name in vertices)
#     len_V = len(vertices) + 1
#
#     weighted_digraph = []
#     input_E = []
#
#     for _, row in df.iterrows():
#         source = row['From']
#         target = row['To']
#         u = V[source]
#         v = V[target]
#         weight = float(row['Effect']) + weight_shift
#         weighted_digraph.append([u, v, weight])
#         input_E.append(f"({u},{v},{weight})")
#
#     input_data = input_V + '\n' + '\n'.join(input_E) + '\n#\n4\ny\n\n\n'
#     print("生成的 GLMY 输入数据：")
#     print(input_data)
#
#     return input_data, weighted_digraph, len_V


def generation_unfiltered_data(csv_path, weight_shift=100, edge_filter='all'):
    df = pd.read_csv(csv_path)

    # 按 edge_type 过滤
    if edge_filter == 'positive':
        df = df[df['edge_type'] == 1]
    elif edge_filter == 'negative':
        df = df[df['edge_type'] == 2]
    elif edge_filter == 'all':
        pass
    else:
        raise ValueError("edge_filter must be one of: 'positive', 'negative', 'all'")

    vertices = sorted(set(df['From']).union(set(df['To'])))
    V = {name: idx + 1 for idx, name in enumerate(vertices)}

    input_V = ','.join(str(V[name]) for name in vertices)
    len_V = len(vertices) + 1

    weighted_digraph = []
    input_E = []

    for _, row in df.iterrows():
        source = row['From']
        target = row['To']
        u = V[source]
        v = V[target]
        weight = float(row['Effect']) + weight_shift
        weighted_digraph.append([u, v, weight])
        input_E.append(f"({u},{v},{weight})")

    input_data = input_V + '\n' + '\n'.join(input_E) + '\n#\n4\ny\n\n\n'
    print("生成的 GLMY 输入数据：")
    print(input_data)

    return input_data, weighted_digraph, len_V


def get_GLMY_result(inputs, exe_file, work_dir, timeout=30):
    process = subprocess.Popen(
        [exe_file],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding='gbk',
        errors='ignore',
        cwd=work_dir
    )

    try:
        if isinstance(inputs, tuple) and len(inputs) > 0:
            input_str = str(inputs[0])
        else:
            input_str = str(inputs)

        print(f"传递给 GLMY 的输入:\n{input_str}")
        output, error = process.communicate(input=input_str, timeout=timeout)

        print("GLMY 输出：")
        print(output)

    except subprocess.TimeoutExpired:
        process.kill()
        output, error = process.communicate()
        return None, "Process timed out.", -1 #超时标记

    return output, error, process.returncode



# 'positive' / 'negative' /all
# CSV → 生成 GLMY 输入 → 调用 GLMY → 检查 homology.json
def process_one_file(csv_path, config):
    file_name = os.path.splitext(os.path.basename(csv_path))[0]
    print(f"\n===== 开始处理: {file_name} =====")

    edge_filter='negative'
    weight_shift=100

    if os.path.exists(config["homology_json_path"]):
        os.remove(config["homology_json_path"])

    inputs = generation_unfiltered_data(
        csv_path=csv_path,
        weight_shift=weight_shift,
        edge_filter=edge_filter,
    )

    output, error, returncode = get_GLMY_result(
        inputs=inputs,
        exe_file=config["exe_file"],
        work_dir=config["work_dir"],
        timeout=config["timeout"]
    )

    print("GLMY return code:", returncode)

    if error:
        print("stderr:")
        print(error[:1000])

    if not os.path.exists(config["homology_json_path"]):
        print(f"{file_name} 未生成 homology.json，跳过。")
        return

    with open(config["homology_json_path"], 'r', encoding='utf-8') as file:
        homology = json.load(file)

    new_json_path = os.path.join(
        os.path.dirname(config["homology_json_path"]),
        f"{file_name}_{edge_filter}.json"
    )

    with open(new_json_path, 'w', encoding='utf-8') as file:
        json.dump(homology, file, ensure_ascii=False, indent=2)

    print(f"已保存为: {new_json_path}")

    edge_filter = 'positive'
    weight_shift = 0

    if os.path.exists(config["homology_json_path"]):
        os.remove(config["homology_json_path"])

    inputs = generation_unfiltered_data(
        csv_path=csv_path,
        weight_shift=weight_shift,
        edge_filter=edge_filter,
    )

    output, error, returncode = get_GLMY_result(
        inputs=inputs,
        exe_file=config["exe_file"],
        work_dir=config["work_dir"],
        timeout=config["timeout"]
    )

    print("GLMY return code:", returncode)

    if error:
        print("stderr:")
        print(error[:1000])

    if not os.path.exists(config["homology_json_path"]):
        print(f"{file_name} 未生成 homology.json，跳过。")
        return

    with open(config["homology_json_path"], 'r', encoding='utf-8') as file:
        homology = json.load(file)

    new_json_path = os.path.join(
        os.path.dirname(config["homology_json_path"]),
        f"{file_name}_{edge_filter}.json"
    )

    with open(new_json_path, 'w', encoding='utf-8') as file:
        json.dump(homology, file, ensure_ascii=False, indent=2)

    print(f"已保存为: {new_json_path}")

    edge_filter = 'all'
    weight_shift = 100

    if os.path.exists(config["homology_json_path"]):
        os.remove(config["homology_json_path"])

    inputs = generation_unfiltered_data(
        csv_path=csv_path,
        weight_shift=weight_shift,
        edge_filter=edge_filter,
    )

    output, error, returncode = get_GLMY_result(
        inputs=inputs,
        exe_file=config["exe_file"],
        work_dir=config["work_dir"],
        timeout=config["timeout"]
    )

    print("GLMY return code:", returncode)

    if error:
        print("stderr:")
        print(error[:1000])

    if not os.path.exists(config["homology_json_path"]):
        print(f"{file_name} 未生成 homology.json，跳过。")
        return

    with open(config["homology_json_path"], 'r', encoding='utf-8') as file:
        homology = json.load(file)

    new_json_path = os.path.join(
        os.path.dirname(config["homology_json_path"]),
        f"{file_name}_{edge_filter}.json"
    )

    with open(new_json_path, 'w', encoding='utf-8') as file:
        json.dump(homology, file, ensure_ascii=False, indent=2)

    print(f"已保存为: {new_json_path}")

    print(f"{file_name} 处理完成。")


def main(config):
    input_dir = config["input_dir"]

    csv_files = sorted([
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.lower().endswith(".csv")
    ])

    if not csv_files:
        print("该文件夹下没有 csv 文件。")
        return

    print(f"共找到 {len(csv_files)} 个 csv 文件。")

    for csv_path in csv_files:
        process_one_file(csv_path, config)

    print("全部文件处理完成。")


if __name__ == "__main__":
    base_root = "D:/Soil3DNetworks/code/GLMY"

    for ii in range(1, 435):
        current_dir = os.path.join(base_root, str(ii))

        CONFIG = {
            "input_dir": current_dir,
            "exe_file": os.path.join(base_root, "GLMY.exe"),
            "work_dir": current_dir,
            "homology_json_path": os.path.join(current_dir, "homology.json"),
            "save_path": current_dir,
            "timeout": 30
        }

        print(f"正在处理第 {ii} 个文件夹：{current_dir}")
        main(CONFIG)