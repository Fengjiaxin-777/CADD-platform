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

st.set_page_config(page_title="QSAR 模型训练", page_icon="🔬", layout="wide")

# ====================== 修复：全局样式 + 中文不乱码 + 修复排版错乱 ======================
st.markdown("""
<style>
    .stApp {
        background-color: #f4f7fa;
    }
    .card-block {
        background-color: #ffffff;
        padding: 24px 28px;
        border-radius: 10px;
        box-shadow: 0 1px 9px rgba(38,50,72,0.07);
        margin-bottom:22px;
        width: 100%;
    }
    h2,h3,h4{
        color:#2c3644;
        font-weight:550;
        font-family: "Microsoft YaHei", sans-serif;
    }
    p,li{
        color:#495057;
        line-height:1.7;
        font-size:15px;
        font-family: "Microsoft YaHei", sans-serif;
    }
    .step-text{
        background:#eef3f9;
        padding:12px 16px;
        border-left:4px solid #597ba8;
        border-radius:4px;
        margin:8px 0;
        width: 100%;
    }
    /* 修复：表格不挤在右边 */
    .stDataFrame {
        width: 100% !important;
    }
    /* 修复：整体宽度铺满 */
    .main .block-container {
        max-width: 1200px;
        width: 100%;
    }
</style>
""", unsafe_allow_html=True)

# 修复：matplotlib 中文
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ====================== 页面介绍 ======================
st.title("定量构效关系(QSAR)建模及可解释性评估")

with st.container(border=False):
    st.markdown('<div class="card-block">', unsafe_allow_html=True)
    st.subheader("一、功能介绍")
    st.markdown("""
    本模块用于构建定量构效关系模型，通过分子结构特征预测生物活性。
    支持多种机器学习算法训练、模型评估、指标可视化与新数据批量预测。
    """)
    st.markdown('</div>', unsafe_allow_html=True)

with st.container(border=False):
    st.markdown('<div class="card-block">', unsafe_allow_html=True)
    st.subheader("二、操作步骤")
    st.markdown("""
    1. 上传带特征与label的CSV文件
    2. 选择机器学习模型并调整参数
    3. 训练模型并查看评估指标
    4. 上传新数据进行批量预测
    """)
    st.markdown('</div>', unsafe_allow_html=True)

# ====================== 特征识别 ======================
def get_feature_columns(df):
    fixed_features = ["MW", "LogP", "TPSA", "HBD", "HBA", "RotB", "QED"]
    morgan_features = [col for col in df.columns if col.startswith("MorganFP")]
    return list(set([c for c in fixed_features if c in df.columns] + morgan_features))

# ====================== 上传数据 ======================
st.markdown('<div class="card-block">', unsafe_allow_html=True)
st.subheader("三、上传训练数据集")
uploaded_files = st.file_uploader("上传CSV文件", type=["csv"], accept_multiple_files=True)

df_combined = None
X_train = X_test = y_train = y_test = scaler = None

if uploaded_files:
    try:
        df_list = [pd.read_csv(file) for file in uploaded_files]
        df_combined = pd.concat(df_list, ignore_index=True)
        st.success("数据上传并合并成功")

        st.dataframe(df_combined.head(10), use_container_width=True)

        feature_cols = get_feature_columns(df_combined)
        label_col = "label" if "label" in df_combined.columns else None

        if len(feature_cols) == 0:
            st.error("未识别到有效特征")
        elif label_col is None:
            st.error("缺少label列")
        else:
            X = df_combined[feature_cols].fillna(0)
            y = df_combined[label_col]
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            X_train, X_test, y_train, y_test = train_test_split(
                X_scaled, y, test_size=0.2, random_state=42, stratify=y
            )
            st.info(f"训练集 {X_train.shape[0]} | 测试集 {X_test.shape[0]}")

    except Exception as e:
        st.error(f"错误：{str(e)}")
st.markdown('</div>', unsafe_allow_html=True)

# ====================== 模型选择 ======================
st.markdown('<div class="card-block">', unsafe_allow_html=True)
st.subheader("四、选择模型与参数")
model_option = st.selectbox("选择模型", [
    "随机森林", "支持向量机", "逻辑回归"
])

params = {"random_state": 42}
if model_option == "随机森林":
    params["n_estimators"] = st.slider("决策树数量", 50, 300, 100)
    params["max_depth"] = st.slider("最大深度", 3, 15, 8)
elif model_option == "支持向量机":
    params["C"] = st.slider("正则化强度", 0.1, 10.0, 1.0)
    params["kernel"] = "rbf"
    params["probability"] = True
else:
    params["C"] = st.slider("正则化强度", 0.1, 10.0, 1.0)
    params["max_iter"] = 500
st.markdown('</div>', unsafe_allow_html=True)

# ====================== 训练 ======================
st.markdown('<div class="card-block">', unsafe_allow_html=True)
st.subheader("五、模型训练")
train_btn = st.button("开始训练")

if train_btn:
    if df_combined is None:
        st.error("请先上传数据")
    else:
        with st.spinner("训练中"):
            if model_option == "随机森林":
                model = RandomForestClassifier(**params)
            elif model_option == "支持向量机":
                model = SVC(**params)
            else:
                model = LogisticRegression(**params)

            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            y_proba = model.predict_proba(X_test)[:, 1]

        acc = accuracy_score(y_test, y_pred)
        pre = precision_score(y_test, y_pred, zero_division=0)
        rec = recall_score(y_test, y_pred, zero_division=0)
        f1 = f1_score(y_test, y_pred, zero_division=0)
        auc = roc_auc_score(y_test, y_proba)

        # 修复：去掉 index=False
        df_perf = pd.DataFrame({
            "指标": ["准确率", "精确率", "召回率", "F1", "AUC"],
            "数值": [acc, pre, rec, f1, auc]
        })
        st.dataframe(df_perf, use_container_width=True)

        # 混淆矩阵
        cm = confusion_matrix(y_test, y_pred)
        fig, ax = plt.subplots(figsize=(3, 2))
        sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
                    xticklabels=["0", "1"],
                    yticklabels=["0", "1"])
        ax.set_xlabel("预测")
        ax.set_ylabel("真实")
        st.pyplot(fig)

        # 保存模型
        buf = io.BytesIO()
        joblib.dump(model, buf)
        buf.seek(0)
        st.download_button("下载模型", buf, "qsar_model.pkl")

        st.session_state['model'] = model
        st.session_state['scaler'] = scaler
        st.session_state['feats'] = feature_cols
st.markdown('</div>', unsafe_allow_html=True)

# ====================== 预测 ======================
st.markdown('<div class="card-block">', unsafe_allow_html=True)
st.subheader("六、新数据预测")
pred_file = st.file_uploader("上传预测数据", type=["csv"])

if pred_file is not None:
    if 'model' not in st.session_state:
        st.warning("请先训练模型")
    else:
        df_pred = pd.read_csv(pred_file)
        feats = get_feature_columns(df_pred)

        if set(feats) != set(st.session_state['feats']):
            st.error("特征不匹配")
        else:
            X_new = df_pred[feats].fillna(0)
            X_new_s = st.session_state['scaler'].transform(X_new)
            y_out = st.session_state['model'].predict(X_new_s)
            df_pred["预测标签"] = y_out

            st.dataframe(df_pred.head(10), use_container_width=True)
            st.download_button("下载结果", df_pred.to_csv(index=False, encoding="utf-8-sig"), "pred.csv")
st.markdown('</div>', unsafe_allow_html=True)
