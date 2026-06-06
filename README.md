# CADD-platform
本平台是一个基于 Streamlit 开发的轻量级 Web 应用，旨在构建药物计算机辅助设计 (CADD) 的一体化工作流。它为药物研发人员提供了一站式解决从分子数清洗质控、分子对接模拟到三维可视化分析的便捷工具，大幅提升初筛效率。

🎯 系统核心功能
1. 数据智能清洗 (Data Quality Control)
理化性质筛选：内置分子量 (MW)、LogP、TPSA、旋转键数 (RB) 等参数的自定义筛选功能。
结构合规审计：自动剔除无法被 RDKit 解析的错误 SMILES 结构并去重。
标准化导出：支持清洗后的高质量分子库一键导出。
2. 自动化分子对接 (Molecular Docking)
物理模拟引擎：集成 AutoDock Vina 引擎，实现蛋白质与小分子的非共价结合自由能精准预测。
参数可控：交互式 Grid Box 设置，支持灵活定义对接活性口袋范围。
批量计算：支持单分子对接与大规模配体库的自动化虚拟筛选。
3. 可交互三维可视化 (Visualization)
3D 渲染视图：内置 py3Dmol，在浏览器内实现蛋白质-配体复合物的实时三维渲染。
结果分析：直观的结合能排行面板，支持一键下载对接报告与 PDB 结构文件包。
🛠 技术栈
前端/工作流: Streamlit
化学信息学: RDKit
分子对接: AutoDock Vina
三维可视化: py3Dmol
数据分析: Pandas, NumPy
🚀 快速开始
环境依赖
项目运行需要如下关键 Python 库，请确保你的环境中已安装：

pip install -r requirements.txt
运行应用
安装依赖后，在根目录下执行：

streamlit run app.py
注：如运行在服务器端，请确保系统已配置 AutoDock Vina 的二进制执行文件路径。

📂 项目结构概览
CADD-platform/
├── app.py                # 项目入口文件
├── pages/                # 业务逻辑页面 (如：数据质控、对接计算)
├── requirements.txt      # 依赖环境配置
├── assets/               # 静态资源、Logo 等
└── README.md             # 项目说明文档
🤝 开发贡献与建议
本项目旨在优化药物设计的工作流，非常欢迎提交 Issue 或 Pull Request 进行功能扩展。

本项目为科研辅助工具，计算结果仅供参考，请结合专业实验数据进行判断。
