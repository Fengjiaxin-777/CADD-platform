# CADD-platform
本平台是一个基于 Streamlit 开发的轻量级药物计算机辅助设计 (CADD) 一体化工作流 Web 应用。它旨在整合药物研发中的关键计算环节，为药物研发人员提供高效、直观的一站式计算和可视化分析工具。

🎯 系统功能模块
平台内置了完善的 CADD 工作流，共分为以下六大核心模块：

1️⃣ 分子质控与物理参数初筛 (1-分子质控与物理参数初筛.py)
格式审计：智能清洗输入的分子数据（如 SMILES），自动去重并剔除无效结构。
物理属性计算：基于 Lipinski 五规则（分子量 MW、LogP、TPSA、旋转键数、氢键给体/受体受限度等）进行初筛。

2️⃣ 分子特征探索 (2-分子特征探索.py)
化学描述符计算：提取小分子的一维/二维/三维分子描述符与特征。
分子指纹构建：支持生成 Morgan/MACCS 等分子指纹，辅助分子的结构相似度分析与聚类。

3️⃣ 成药性与 ADMET 药代动力学评估 (3-成药性与 ADMET 药代动力学评估.py)
ADMET 预测模型：快速预测候选药物分子的吸收 (Absorption)、分布 (Distribution)、代谢 (Metabolism)、排泄 (Excretion) 和毒性 (Toxicity) 性质。
成药性初筛：识别高毒性或药代性质极差的分子，降低先导化合物的后期开发风险。

4️⃣ 定量构效关系建模及可解释性评估 (4-定量构效关系建模及可解释性评估.py)
QSAR 建模：利用主流机器学习算法（如随机森林、支持向量机等）建立分子结构与生物活性的定量构效模型。
可解释性分析：评估分子特征重要性，揭示影响活性的关键子结构或化学物理参数。

5️⃣ 三维配体受体物理对接 (5-分子对接.py)
分子对接解算：内部集成 AutoDock Vina 引擎进行高精度的受体-配体结合亲和力模拟。
位点网格盒设置：交互式确定活性口袋 Grid Box 的三维重心坐标及空间尺寸。
3D 交互式可视化：利用 py3Dmol 与 stmol 在网页中生成高质量的受体与配体结合复合物的 3D 可视化交互图。

6️⃣ 数据解析与可视化 (6-数据解析与可视化.py)
多维数据绘图：内置 Matplotlib、Seaborn 和 Plotly 动态分析图表，对筛选数据和对接得分进行深度统计学解析。
打包导出：对接排行结果支持一键导出为 CSV 报告，且包含所有优化构象的 PDB 压缩包下载。

🛠 技术栈与依赖
开发框架: Streamlit
化学信息学核心: RDKit
分子对接解算: AutoDock Vina
三维分子渲染: py3Dmol, stmol
科学计算 & 机器学习: Pandas, NumPy, Scikit-learn, joblib
数据可视化: Matplotlib, Seaborn, Plotly

🚀 快速部署与本地运行
1. 克隆代码仓库到本地
git clone https://github.com/Fengjiaxin-777/CADD-platform.git
cd CADD-platform
2. 安装系统底层依赖
若在 Linux 或 Streamlit Cloud 环境部署，本平台可能需要 openbabel 等依赖以转换不同格式的分子文件。请确保系统包含 packages.txt 中的依赖项（例如 openbabel）。

3. 安装 Python 依赖包
pip install -r requirements.txt
4. 运行 Streamlit 服务器
streamlit run Home.py
启动后在浏览器中打开：http://localhost:8501。

📂 仓库项目结构
CADD-platform/

├── Home.py                  # 平台主页 (Portal)

├── requirements.txt         # Python 依赖清单

├── packages.txt             # Linux 系统底层包依赖清单 (如 openbabel)

├── data/                    # 示例数据与缓存文件夹

└── pages/                   # 系统子功能控制页

    ├── 1-分子质控与物理参数初筛.py
    
    ├── 2-分子特征探索.py
    
    ├── 3-成药性与 ADMET 药代动力学评估.py
    
    ├── 4-定量构效关系建模及可解释性评估.py
    
    ├── 5-分子对接.py
    
    └── 6-数据解析与可视化.py
    
✍️ 说明
本平台中的预测结果（如 ADMET/QSAR）为机器学习模拟值，不能完全代替体内/体外湿实验验证。

3D 可视化时，如果找不到分子结构，请确保对接过程已顺利生成对应的 PDB 文件并且未被工作区清理。

