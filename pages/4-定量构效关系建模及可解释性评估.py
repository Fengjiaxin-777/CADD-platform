# -*- coding: utf-8 -*-
import streamlit as st
import joblib
import io
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import (accuracy_score, precision_score, recall_score,
                             f1_score, roc_auc_score, confusion_matrix)
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')

# 解决特定环境下可能导致的输出流编码问题
import sys
if sys.stdout.encoding != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# 字体配置，多系统兼容方案（特别针对 Linux/Streamlit Cloud 部署环境加入了 WenQuanYi）
def set_matplot_zh_font():
    plt.rcParams['font.sans-serif'] = [
        'SimHei',             # Windows 黑体
        'Microsoft YaHei',    # Windows 微软雅黑
        'Heiti TC',           # macOS 繁体黑体
        'Arial Unicode MS',   # macOS 兼容中文
        'WenQuanYi Micro Hei',# Linux 下的开源中文字体
        'DejaVu Sans',        # 备用英文
        'sans-serif'
    ]
    plt.rcParams['axes.unicode_minus'] = False # 正常显示负号

set_matplot_zh_font()

st.set_page_config(
    page_title="QSAR 模型训练工具",
    layout="wide"
)

st.markdown("""
<style>
    .main { background-color: #f8fafc; }
    .blue-banner {
        background: linear-gradient(135deg, #0f172a 0%, #1e3a8a 100%);
        color: #ffffff;
        padding: 30px;
        border-radius: 12px;
        margin-bottom: 25px;
    }
    .blue-banner h1 { color: #ffffff !important; margin: 0 0 8px 0 !important; font-size: 28px; }
    .blue-banner p { color: #cbd5e1 !important; margin: 0 !important; font-size: 14px; }
    
    /* 统一的区块标题风格（浅灰底色+深蓝左色条） */
    .section-header {
        background-color: #f1f5f9;
        padding: 10px 16px;
        border-radius: 8px;
        color: #0f172a;
        font-weight: 700;
        font-size: 16px;
        margin-bottom: 15px;
        border-left: 5px solid #1e3a8a;
    }
    .report-card {
        background-color: #ffffff;
        border-left: 4px solid #1e3a8a;
        padding: 20px;
        border-radius: 8px;
        border: 1px solid #e2e8f0;
        margin-top: 15px;
        margin-bottom: 15px;
    }
    
    /* 强制重置 Streamlit 原生 border 容器的边框和内边距，统一风格 */
    div[data-testid="stVerticalBlockBorderLine"] {
        background-color: #ffffff !important;
        border: 1px solid #e2e8f0 !important;
        border-radius: 12px !important;
        padding: 22px !important;
        box-shadow: 0 1px 3px 0 rgba(0, 0, 0, 0.05) !important;
        margin-bottom: 20px !important;
    }
    .stButton>button {
        background-color: #1e3a8a !important;
        color: white !important;
        border-radius: 8px !important;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="blue-banner">
    <h1>QSAR 机器学习模型训练平台</h1>
    <p>支持多分子特征表格合并，快速构建并评测经典 QSAR 分类模型，实现药物分子的活性定量分析与筛选。</p>
</div>
""", unsafe_allow_html=True)

# 初始化 session 状态以维持持久化变量
if "X_train" not in st.session_state:
    st.session_state["X_train"] = None
    st.session_state["X_test"] = None
    st.session_state["y_train"] = None
    st.session_state["y_test"] = None
    st.session_state["scaler"] = None
    st.session_state["feature_cols"] = None
    st.session_state["X_raw_mean"] = None
    st.session_state["trained_model"] = None

# ==========================================
# 布局配置 1：格式要求与数据集上传并排
# ==========================================
with st.container(border=True):
    st.markdown('<div class="section-header">数据导入与格式校验</div>', unsafe_allow_html=True)
    col_help, col_upload = st.columns([0.45, 0.55])
    
    with col_help:
        st.markdown("""
        <div style="background-color: #f8fafc; border: 1px dashed #cbd5e1; padding: 12px; border-radius: 8px; font-size: 13px; color: #475569; line-height: 1.6;">
            <b>特征列识别规则</b>：<br>
            1. 默认识别以下理化特征列：<code>Lipinski_Violations, QED, MW, LogP, TPSA, HBD, HBA, RotatableBonds, RingCount, AromaticRingCount, HeteroatomCount, FractionCSP3, MaxPartialCharge, MinPartialCharge, Polarizability</code><br>
            2. 指纹前缀匹配：所有以 <code>MorganFP_</code> 开头的字段列均将作为指纹特征纳入计算。<br>
            <b>分类标签列</b>：<br>
            必须存在一列名为 <code>label</code> 的数字标签，值仅限为 <code>0</code>（阴性）或 <code>1</code>（阳性）。
        </div>
        """, unsafe_allow_html=True)
        
    with col_upload:
        uploaded_files = st.file_uploader("选择一个或多个 CSV 文件进行训练：", type=["csv"], accept_multiple_files=True)

# 特征提取辅助函数
def get_feature_columns(df):
    fixed_features = [
        "Lipinski_Violations", "QED", "MW", "LogP", "TPSA", "HBD", "HBA",
        "RotatableBonds", "RingCount", "AromaticRingCount", "HeteroatomCount",
        "FractionCSP3", "MaxPartialCharge", "MinPartialCharge", "Polarizability",
        "LipinskiViolations"
    ]
    morgan_features = [col for col in df.columns if col.startswith("MorganFP")]
    valid_features = [col for col in fixed_features if col in df.columns] + morgan_features
    return list(set(valid_features))

# 校验并提取载入数据
df_combined = None
if uploaded_files:
    try:
        df_list = []
        for file in uploaded_files:
            df_temp = pd.read_csv(file)
            df_list.append(df_temp)
            
        df_combined = pd.concat(df_list, ignore_index=True)
        
        with st.container(border=True):
            st.markdown('<div class="section-header">数据集校对与预览</div>', unsafe_allow_html=True)
            st.success(f"成功导入并自动合并 {len(uploaded_files)} 个数据集文件。")
            
            with st.expander("点击展开合并后原始数据集预览"):
                st.dataframe(df_combined.head(8), use_container_width=True)
                st.caption(f"当前合并表尺寸: {df_combined.shape[0]} 行 * {df_combined.shape[1]} 列")
            
            feature_cols = get_feature_columns(df_combined)
            label_col = "label" if "label" in df_combined.columns else None

            if len(feature_cols) == 0:
                st.error("数据校验未通过：未在数据表内识别到符合标准的理化特征或分子特征指纹列。")
                df_combined = None
            elif label_col is None:
                st.error("数据校验未通过：未找到指定的样本活性标签列 label。")
                df_combined = None
            elif not df_combined[label_col].isin([0, 1]).all():
                st.error("数据校验未通过：标签列 label 中的元素仅限于分类数字 0 与 1。")
                df_combined = None
            else:
                st.info(f"数据结构自动识别完成：已筛选出可用于 QSAR 建模的特征列共 {len(feature_cols)} 个。")
                
                # 分割特征与标签并持久化处理
                X = df_combined[feature_cols]
                y = df_combined[label_col]
                
                X_filled = X.fillna(X.mean())
                scaler = StandardScaler()
                X_scaled = scaler.fit_transform(X_filled)
                
                X_train, X_test, y_train, y_test = train_test_split(
                    X_scaled, y, test_size=0.2, random_state=42, stratify=y
                )
                
                # 存入 Session State 中，确保按钮操作不触发重新上传和复位
                st.session_state["X_train"] = X_train
                st.session_state["X_test"] = X_test
                st.session_state["y_train"] = y_train
                st.session_state["y_test"] = y_test
                st.session_state["scaler"] = scaler
                st.session_state["feature_cols"] = feature_cols
                st.session_state["X_raw_mean"] = X.mean()
                
                st.markdown(f"""
                <div style="background-color: #f0fdf4; border-left: 4px solid #16a34a; padding: 12px; border-radius: 6px; font-size: 13px; color: #15803d;">
                    数据集分割完成：用以训练的样本数为 <b>{X_train.shape[0]}</b> 个，测试与验证的特征样本数为 <b>{X_test.shape[0]}</b> 个。
                </div>
                """, unsafe_allow_html=True)
            
    except Exception as e:
        st.error(f"数据集解析或自动关联异常，原因: {str(e)}")
        df_combined = None

# ==========================================
# 布局配置 2：模型选择、调参与训练对比
# ==========================================
if st.session_state["X_train"] is not None:
    with st.container(border=True):
        st.markdown('<div class="section-header">QSAR 分类模型算法参数调控</div>', unsafe_allow_html=True)
        
        col_sel, col_param = st.columns([0.4, 0.6])
        
        with col_sel:
            model_option = st.selectbox(
                "选择算法类型：",
                ("随机森林 (Random Forest)", "支持向量机 (SVM)", "逻辑回归 (Logistic Regression)")
            )
        
        with col_param:
            params = {}
            if model_option == "随机森林 (Random Forest)":
                params["n_estimators"] = st.slider("树的个数 (n_estimators)", 50, 500, 100, 10)
                params["max_depth"] = st.slider("最大树深 (max_depth)", 3, 20, 10)
                params["random_state"] = 42
            elif model_option == "支持向量机 (SVM)":
                params["C"] = st.slider("惩罚权重参数 C", 0.01, 10.0, 1.0, 0.1)
                params["kernel"] = st.radio("核函数类型：", ["linear", "rbf"], index=0, horizontal=True)
                params["random_state"] = 42
                params["probability"] = True
            elif model_option == "逻辑回归 (Logistic Regression)":
                params["C"] = st.slider("反正则化强度参数 C", 0.01, 10.0, 1.0, 0.1)
                params["max_iter"] = st.slider("求解器寻优最大迭代次数", 100, 1000, 200, 100)
                params["random_state"] = 42

        st.markdown("<br>", unsafe_allow_html=True)
        train_btn = st.button("运行当前配置单模型训练", type="primary", use_container_width=True)
        
        if train_btn:
            with st.spinner("当前分类算法参数评估学习中，请稍候..."):
                X_train = st.session_state["X_train"]
                X_test = st.session_state["X_test"]
                y_train = st.session_state["y_train"]
                y_test = st.session_state["y_test"]
                
                if model_option == "随机森林 (Random Forest)":
                    model = RandomForestClassifier(**params)
                elif model_option == "支持向量机 (SVM)":
                    model = SVC(**params)
                elif model_option == "逻辑回归 (Logistic Regression)":
                    model = LogisticRegression(**params)
                    
                model.fit(X_train, y_train)
                st.session_state["trained_model"] = model
                
                y_pred = model.predict(X_test)
                y_pred_proba = model.predict_proba(X_test)[:, 1]
                
                # 计算各项测试统计指标
                accuracy = accuracy_score(y_test, y_pred)
                precision = precision_score(y_test, y_pred)
                recall = recall_score(y_test, y_pred)
                f1 = f1_score(y_test, y_pred)
                auc = roc_auc_score(y_test, y_pred_proba)
                
                st.success("模型学习与验证集分类判定工作全部结束。")
                
                # 报告展示区
                meta_col1, meta_col2, meta_col3, meta_col4, meta_col5 = st.columns(5)
                meta_col1.metric("分类准确度 (Accuracy)", f"{accuracy:.3f}")
                meta_col2.metric("精确度 (Precision)", f"{precision:.3f}")
                meta_col3.metric("灵敏度/召回率 (Recall)", f"{recall:.3f}")
                meta_col4.metric("F1 综合分数", f"{f1:.3f}")
                meta_col5.metric("受试者特征面积 (AUC)", f"{auc:.3f}")
                
                # 混淆矩阵热图
                st.markdown("<br>", unsafe_allow_html=True)
                st.markdown("**测试集混淆分类矩阵分布：**")
                cm = confusion_matrix(y_test, y_pred)
                fig_cm, ax_cm = plt.subplots(figsize=(4.5, 3.5))
                sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
                            xticklabels=["阴性(Negative/0)", "阳性(Positive/1)"],
                            yticklabels=["阴性(Negative/0)", "阳性(Positive/1)"], ax=ax_cm)
                ax_cm.set_ylabel("实际标记类别 (Real Label)")
                ax_cm.set_xlabel("预估判定标签 (Predicted Label)")
                plt.tight_layout()
                st.pyplot(fig_cm)
                plt.close(fig_cm)
                
                # 模型持久化导出下载
                model_buf = io.BytesIO()
                joblib.dump(model, model_buf)
                model_buf.seek(0)
                st.download_button(
                    label="导出本轮构建的分类模型权重二进制包 (.pkl)",
                    data=model_buf,
                    file_name="qsar_classifier_model.pkl",
                    mime="application/octet-stream",
                    use_container_width=True
                )

    # 2.2 横向评估模型性能对比板块
    with st.container(border=True):
        st.markdown('<div class="section-header">QSAR 算法模型横向评测与选型对比</div>', unsafe_allow_html=True)
        st.caption("注：横向对比时，系统将会统一调用基准通用参数在后台对三种不同的算法进行同时训练，并分析对比结果。")
        
        run_comparison = st.button("启动主流模型分类效能横向对比评估", use_container_width=True)
        
        if run_comparison:
            with st.spinner("横向预测效能对比训练中..."):
                X_train = st.session_state["X_train"]
                X_test = st.session_state["X_test"]
                y_train = st.session_state["y_train"]
                y_test = st.session_state["y_test"]
                
                comparative_classifiers = {
                    "随机森林 (Random Forest)": RandomForestClassifier(n_estimators=150, max_depth=10, random_state=42),
                    "支持向量机 (SVM)": SVC(C=1.0, kernel='rbf', probability=True, random_state=42),
                    "逻辑回归 (Logistic Regression)": LogisticRegression(C=1.0, max_iter=300, random_state=42)
                }
                
                results_records = []
                for name, clf in comparative_classifiers.items():
                    clf.fit(X_train, y_train)
                    y_val_pred = clf.predict(X_test)
                    y_val_proba = clf.predict_proba(X_test)[:, 1]
                    
                    results_records.append({
                        "机器学习算法": name,
                        "准确率 (Accuracy)": round(accuracy_score(y_test, y_val_pred), 3),
                        "精确率 (Precision)": round(precision_score(y_test, y_val_pred), 3),
                        "召回率 (Recall)": round(recall_score(y_test, y_val_pred), 3),
                        "F1得分 (F1-score)": round(f1_score(y_test, y_val_pred), 3),
                        "AUC曲线下面积 (ROC-AUC)": round(roc_auc_score(y_test, y_val_proba), 3)
                    })
                    
                compare_df = pd.DataFrame(results_records)
                st.markdown("**三种经典 QSAR 分类算法评估指标明细：**")
                st.dataframe(compare_df, use_container_width=True)
                
                st.markdown("<br>", unsafe_allow_html=True)
                st.markdown("**不同模型指标效果直观对比图表：**")
                
                plot_df = compare_df.melt(id_vars="机器学习算法", var_name="评价指标", value_name="评分")
                
                fig_comp, ax_comp = plt.subplots(figsize=(10, 4.5))
                sns.barplot(data=plot_df, x="评价指标", y="评分", hue="机器学习算法", palette="Blues_r", ax=ax_comp)
                ax_comp.set_ylim(0, 1.1)
                ax_comp.set_xlabel("不同维度评定指标特征", fontsize=10)
                ax_comp.set_ylabel("最终判定评分", fontsize=10)
                ax_comp.spines['top'].set_visible(False)
                ax_comp.spines['right'].set_visible(False)
                plt.legend(frameon=True, loc='lower right')
                plt.tight_layout()
                st.pyplot(fig_comp)
                plt.close(fig_comp)

    # 2.3 上传新无标签数据集批量预测
    with st.container(border=True):
        st.markdown('<div class="section-header">未知活性的新分子配体预测判定</div>', unsafe_allow_html=True)
        pred_file = st.file_uploader("导入待分类判定并无标签的数据文件（SMILES或具有同等理化常数特征列的 CSV）：", type=["csv"])
        
        if pred_file is not None:
            if st.session_state["trained_model"] is None:
                st.warning("提示：请先在上方确定分类模型参数，点击“运行当前配置单模型训练”完成模型准备后，方可启动本分类器。")
            else:
                df_pred_raw = pd.read_csv(pred_file)
                pred_feats = get_feature_columns(df_pred_raw)
                feature_cols = st.session_state["feature_cols"]
                
                if set(pred_feats) != set(feature_cols):
                    st.error(f"字段指纹及描述特征列不一致！模型所需特征数：{len(feature_cols)}。当前诊断出特征数：{len(pred_feats)}。请重新比对表头名称。")
                else:
                    with st.spinner("系统导入预处理性质，并计算最终决策类别..."):
                        mean_val = st.session_state["X_raw_mean"]
                        X_pred = df_pred_raw[pred_feats].fillna(mean_val)
                        
                        scaler = st.session_state["scaler"]
                        model = st.session_state["trained_model"]
                        
                        X_pred_scaled = scaler.transform(X_pred)
                        y_pred_out = model.predict(X_pred_scaled)
                        
                        df_out_res = df_pred_raw.copy()
                        df_out_res["分类预测值 (Predicted_Label)"] = y_pred_out
                        
                        st.success("批量预测诊断完成。")
                        st.markdown("**预测分类得到的化合物属性列表前缀记录（前10行记录）：**")
                        st.dataframe(df_out_res.head(10), use_container_width=True)
                        
                        res_csv = df_out_res.to_csv(index=False).encode("utf-8-sig")
                        st.markdown("<br>", unsafe_allow_html=True)
                        st.download_button(
                            label="导出分子库新样本筛选判断类别表格结果 (CSV)",
                            data=res_csv,
                            file_name="qsar_discovery_prediction_report.csv",
                            mime="text/csv",
                            use_container_width=True
                        )
else:
    st.info("提示：请先在页面顶端导入经过预处理合并的配体 CSV 数据集进行模型初始化。")
