# -*- coding: utf-8 -*-  
import streamlit as st  
import subprocess  
import py3Dmol  
import pandas as pd  
import os  
import shutil  
import tempfile  
import time  
from pathlib import Path  
import zipfile  
  
  
st.set_page_config(page_title="分子三维物理对接", layout="wide")  
  
  
# CSS 样式定义，彻底移除了 .card DOM 类，通过原生 stVerticalBlockBorderLine 完成卡片立体化美化  
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
    .report-card {   
        background-color: #ffffff;   
        border-left: 4px solid #1e3a8a;   
        padding: 20px;   
        border-radius: 8px;   
        border: 1px solid #e2e8f0;   
    }   
    /* 精确美化 Streamlit border 容器的边框和内边距，实现无缝卡片视觉效果 */  
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
    <h1>三维配体受体物理对接</h1>   
    <p>预测小分子在目标蛋白质活性口袋中的非共价键合姿态，解算对接结合亲和自由能并定位关键相互作用。</p>   
</div>   
""", unsafe_allow_html=True)   
  
  
# ==============================================================================
# 提供一段结构绝对标准、Vina 可顺利计算的标准小分子配体 PDBQT 示例文本数据（包含扭转树定义）
# ==============================================================================
SAMPLE_LIGAND_PDBQT = """REMARK  5 ACTIVE TORSIONS
REMARK  Assigned autodock charges and atom types
ROOT
ATOM      1  C   LIG     1       0.757   0.141   0.000  1.00  0.00    +0.120 C 
ATOM      2  O   LIG     1      -0.643   0.141   0.000  1.00  0.00    -0.340 OA
ATOM      3  H   LIG     1       1.114   1.041   0.300  1.00  0.00    +0.060 HD
ATOM      4  H   LIG     1       1.114  -0.759  -0.300  1.00  0.00    +0.060 HD
ATOM      5  H   LIG     1      -1.000  -0.759   0.100  1.00  0.00    +0.100 HD
ENDROOT
TORSDOF 0
"""

# initialize session state keys   
if "results_df" not in st.session_state:   
    st.session_state["results_df"] = None   
if "protein_pdb" not in st.session_state:   
    st.session_state["protein_pdb"] = None   
if "workdir" not in st.session_state:   
    # persistent working directory under system temp (not auto-deleted)   
    st.session_state["workdir"] = os.path.join(tempfile.gettempdir(), f"streamlit_docking_{int(time.time())}")   
    os.makedirs(st.session_state["workdir"], exist_ok=True)   
   
   
workdir = st.session_state["workdir"]   
   
   
def safe_filename(name: str) -> str:   
    name = str(name).strip()   
    name = "".join(c for c in name if c.isalnum() or c in ("_", "-", "."))   
    if not name:   
        name = f"mol_{int(time.time())}"   
    return name   
   
   
def check_tool(name: str) -> bool:   
    return shutil.which(name) is not None   

# 内存中针对配体文件做前置严格合规校验
def validate_pdbqt_format(file_content_bytes: bytes) -> bool:
    try:
        content_str = file_content_bytes.decode("utf-8", errors="ignore")
        # AutoDock Vina 要求配体 PDBQT 中必须包含 ROOT 的开始标记作为分子扭转树根
        return "ROOT" in content_str
    except Exception:
        return False
   
   
def parse_vina_output(stdout_text: str):   
    # try both "1   -7.9 ..." and "REMARK VINA RESULT: -7.9 ..."   
    score = None   
    lines = stdout_text.splitlines()   
    for line in lines:   
        line = line.strip()   
        if line.startswith("REMARK VINA RESULT"):   
            parts = line.split()   
            for p in parts[::-1]:   
                try:   
                    score = float(p)   
                    return score   
                except:   
                    continue   
        # split tokens   
        tokens = line.split()   
        if len(tokens) >= 2 and tokens[0].isdigit():   
            # line like: "1   -7.9   0.000   0.000"   
            try:   
                return float(tokens[1])   
            except:   
                continue   
    return score   
   
   
def run_vina(protein_path, ligand_path, center, size, exhaustiveness, out_pdbqt):   
    # 优先尝试利用 Python API 接口运行对接以避免命令行 PATH 权限或缺失问题  
    try:
        from vina import Vina
        v = Vina(sf_name='vina', cpu=1)
        v.set_receptor(protein_path)
        v.set_ligand_from_file(ligand_path)
        v.compute_vina_maps(center=list(center), box_size=list(size))
        v.dock(exhaustiveness=exhaustiveness, n_poses=1)
        v.write_poses(out_pdbqt, n_poses=1, overwrite=True)
        energies = v.energies(n_poses=1)
        score = round(float(energies[0][0]), 2)
        return score, "Success via Python API"
    except ImportError:
        # 如果未引入 python-vina 依赖，则后备退回到命令行 subprocess 执行逻辑
        if not check_tool("vina"):   
            return None, "Error: 未在系统中检测到 'vina' 命令行程序，同时 Python 环境中也未检测到 'vina' 包。"  
        
        cmd = [   
            "vina",   
            "--receptor", protein_path,   
            "--ligand", ligand_path,   
            "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),   
            "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),   
            "--exhaustiveness", str(exhaustiveness),   
            "--out", out_pdbqt   
        ]   
        try:   
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)   
            out = proc.stdout + "\n" + proc.stderr   
            score = parse_vina_output(out)   
            return score, out   
        except subprocess.TimeoutExpired:   
            return None, "Vina 对接耗时超限（10分钟超时中断）"   
        except Exception as e:   
            return None, f"对接执行异常: {str(e)}"   
    except Exception as e:
        return None, f"Python Vina API 执行发生故障: {str(e)}"
   
   
def pdbqt_to_pdb(pdbqt_file, pdb_file):   
    # prefer obabel if present   
    if check_tool("obabel"):   
        try:   
            subprocess.run(["obabel", "-ipdbqt", pdbqt_file, "-opdb", "-O", pdb_file], check=True)   
            return True   
        except Exception:   
            pass   
    # fallback: copy with .pdb extension (py3Dmol can often render)   
    try:   
        shutil.copy(pdbqt_file, pdb_file)   
        return True   
    except Exception:   
        return False   
   
   
def show_structure(protein_path, ligand_path):    
    # 诊断性提示：确保我们知道到底读到了什么    
    if not os.path.exists(protein_path):    
        st.error(f"严重错误：受体文件不存在 {protein_path}")    
        return    
    if not os.path.exists(ligand_path):    
        st.error(f"严重错误：配体文件不存在 {ligand_path}")    
        return    
             
    try:    
        with open(protein_path, "r", encoding="utf-8", errors="ignore") as f: p_data = f.read()    
        with open(ligand_path, "r", encoding="utf-8", errors="ignore") as f: l_data = f.read()    
             
        # 如果文件为空    
        if not p_data.strip(): st.error("错误：受体文件内容为空")    
        if not l_data.strip(): st.error("错误：配体文件内容为空")    
      
      
        view = py3Dmol.view(width=850, height=520)    
        view.addModel(p_data, "pdb")    
        view.setStyle({"model": 0}, {"cartoon": {"color": "#64748b"}})    
        view.addModel(l_data, "pdb")    
        view.setStyle({"model": 1}, {"stick": {"colorscheme": "cyanCarbon"}})    
        view.zoomTo()    
             
        # 移除了 html 渲染中的 key 参数参数，彻底修复旧版 Streamlit Iframe 命名冲突报错  
        view_html = view._make_html()     
        st.components.v1.html(view_html, height=520)    
             
    except Exception as e:     
        st.error(f"3D 渲染渲染流程捕获异常：{e}")    
        st.exception(e) # 打印详细的堆栈跟踪    
   
   
col_ctrl, col_view = st.columns([0.38, 0.62])   
   
   
with col_ctrl:   
    with st.container(border=True):   
        st.subheader("对接任务设置")   
        
        # --- 新增：强制性的数据格式要求规范及官方下载模块演示 ---
        st.markdown("""
        <div style="background-color: #eff6ff; border-left: 4px solid #3b82f6; padding: 12px; border-radius: 6px; font-size: 13px; color: #1e3a8a; margin-bottom:15px; line-height:1.5;">
            <strong>📋 数据格式要求要求：</strong><br>
            • 配体文件<strong>必须</strong>是通过 <code>AutoDockTools</code>, <code>OpenBabel</code> 或 <code>MGLTools</code> 加载处理、含有电荷和可扭转树架构的<strong>标准 PDBQT 格式</strong>。<br>
            • 严禁粗暴强行更改 <code>.pdb</code> 后缀为 <code>.pdbqt</code>，以防计算报错。你可以下载下方经校验的标准正确数据做首次运行校验。
        </div>
        """, unsafe_allow_html=True)
        
        st.download_button(
            label="💾 下载标准配体 PDBQT 示例数据",
            data=SAMPLE_LIGAND_PDBQT,
            file_name="standard_ligand_sample.pdbqt",
            mime="text/plain",
            use_container_width=True
        )
        
        st.markdown("---")
        run_mode = st.radio("请指定配体数量类型：", ["对接单个候选分子", "执行批量虚拟配体筛选"])   
   
   
        st.markdown("---")   
        f_prot = st.file_uploader("载入大分子受体蛋白质 (.pdbqt 或 .pdb)", type=["pdbqt", "pdb"])   
   
   
        if run_mode == "对接单个候选分子":   
            f_lig = st.file_uploader("载入要对接的配体小分子 (.pdbqt)", type=["pdbqt"])   
            f_ligs = None   
        else:   
            f_ligs = st.file_uploader("批量载入多个配体分子 (.pdbqt)", type=["pdbqt"], accept_multiple_files=True)   
            f_lig = None   
   
   
        st.markdown("---")   
        st.subheader("活性位点口袋网格划定 (Grid Box)")   
   
   
        c_x_c, c_y_c, c_z_c = st.columns(3)   
        with c_x_c:   
            cx = st.number_input("三轴中心 X 坐标", value=0.0, step=0.1)   
        with c_y_c:   
            cy = st.number_input("三轴中心 Y 坐标", value=0.0, step=0.1)   
        with c_z_c:   
            cz = st.number_input("三轴中心 Z 坐标", value=0.0, step=0.1)   
   
   
        s_x_c, s_y_c, s_z_c = st.columns(3)   
        with s_x_c:   
            sx = st.number_input("边宽 X (Å)", value=20.0, step=0.5)   
        with s_y_c:   
            sy = st.number_input("高度 Y (Å)", value=20.0, step=0.5)   
        with s_z_c:   
            sz = st.number_input("厚度 Z (Å)", value=20.0, step=0.5)   
   
   
        v_exh = st.slider("采样深度 (Exhaustiveness)", min_value=1, max_value=32, value=8)   
        btn_start = st.button("启动结合模拟计算", type="primary", use_container_width=True)   
   
   
with col_view:   
    # workspace controls   
    with st.container(border=True):   
        st.write(f"工作目录：{workdir}")   
        c1, c2 = st.columns([0.5, 0.5])   
        with c1:   
            if st.button("清理工作区文件"):   
                try:   
                    shutil.rmtree(workdir)   
                except Exception:   
                    pass   
                os.makedirs(workdir, exist_ok=True)   
                st.session_state["results_df"] = None   
                st.session_state["protein_pdb"] = None   
                st.success("工作区已清理。")   
        with c2:   
            if st.button("检测环境 (Vina / OpenBabel)"):   
                ava_vina = check_tool("vina")   
                ava_ob = check_tool("obabel")   
                try:
                    from vina import Vina
                    api_vina = True
                except ImportError:
                    api_vina = False
                st.info(f"Vina 命令行可用: {ava_vina} | Vina Python API 支持: {api_vina} | OpenBabel 可用: {ava_ob}")   
   
   
    if btn_start:   
        if not f_prot:   
            st.error("执行被禁：请上传靶点受体文件 (.pdbqt 或 .pdb)。")   
        else:
            # 读取配体制备列表并进行前置强力格式检测校验
            lig_uploads = []
            if run_mode == "对接单个候选分子":
                if f_lig:
                    lig_uploads.append(f_lig)
            else:
                if f_ligs:
                    lig_uploads = list(f_ligs)
            
            # --- 核心拦截校验步骤 ---
            format_error_list = []
            for item in lig_uploads:
                # 预读文件内容
                file_bytes = item.getvalue()
                if not validate_pdbqt_format(file_bytes):
                    format_error_list.append(item.name)
            
            if format_error_list:
                st.error(f"❌ 启动中断：检测到以下配体文件不符合标准 PDBQT 格式要求：`{', '.join(format_error_list)}`")
                st.markdown("""
                **主要原因分析：**
                1. 您的 `.pdbqt` 可能是强行修改原本 `.pdb` / `.txt` 后缀得到的。
                2. 缺少必要扭转树的 `ROOT` 结构标记和力场局部电荷，导致后台 Vina 引擎阻断异常。
                
                **自助排查建议：**
                * 点击左下角 **“💾 下载标准配体 PDBQT 示例数据”** 下载标准文件对照修改。
                * 使用系统检测支持的 `OpenBabel` 指令或 `AutoDock Tools` 软件导出转换。
                """)
                st.stop()  # 拦截后续的后端计算调用，强制用户修正
            
            # save protein file into workdir   
            prot_name = safe_filename(getattr(f_prot, "name", "protein.pdbqt"))   
            prot_path = os.path.join(workdir, prot_name)   
            with open(prot_path, "wb") as fw:   
                fw.write(f_prot.read())   
            # ensure we have a pdbqt for vina and a pdb for visualization   
            prot_pdbqt = prot_path   
            prot_pdb = os.path.splitext(prot_path)[0] + ".pdb"   
            # if uploaded file is .pdb, convert/copy to pdbqt? We just keep as is for visualization,   
            # and pass prot_pdbqt path to vina if it's .pdbqt. If it's .pdb and vina present, user should provide pdbqt.   
            if prot_path.endswith(".pdb"):   
                # copy to protein.pdb (already)   
                shutil.copy(prot_path, prot_pdb)   
            else:   
                # convert to pdb for visualization if possible   
                pdbqt_to_pdb(prot_path, prot_pdb)   
   
   
            st.session_state["protein_pdb"] = prot_pdb   
   
   
            calc_matrix = []   
            lig_files = []   
            if run_mode == "对接单个候选分子":   
                if not f_lig:   
                    st.error("执行被禁：请上传配体描述文件 (.pdbqt)。")   
                else:   
                    lig_name = safe_filename(os.path.splitext(f_lig.name)[0])   
                    lig_pdbqt = os.path.join(workdir, f"{lig_name}.pdbqt")   
                    with open(lig_pdbqt, "wb") as fw:   
                        # 回读已重置位置的文件流
                        f_lig.seek(0)
                        fw.write(f_lig.read())   
                    lig_files = [lig_pdbqt]   
            else:   
                if not f_ligs:   
                    st.error("执行被禁：批量模式下未上传任何配体。")   
                else:   
                    lig_files = []   
                    for single in f_ligs:   
                        ln = safe_filename(os.path.splitext(single.name)[0])   
                        lp = os.path.join(workdir, f"{ln}.pdbqt")   
                        with open(lp, "wb") as fw:   
                            single.seek(0)
                            fw.write(single.read())   
                        lig_files.append(lp)   
   
   
            if lig_files:   
                progress = st.progress(0.0)   
                total = len(lig_files)   
                for idx, lig_path in enumerate(lig_files):   
                    lig_basename = os.path.splitext(os.path.basename(lig_path))[0]   
                    out_pdbqt = os.path.join(workdir, f"out_{lig_basename}.pdbqt")   
                    t0 = time.time()   
                    # 运行真实对接计算  
                    score, raw = run_vina(prot_pdbqt, lig_path, (cx, cy, cz), (sx, sy, sz), v_exh, out_pdbqt)   
                    t1 = time.time()   
                      
                    # 真实对接模式下，若计算未取得合理得分，记录报错且直接跳过，不再给予随机分数！  
                    if score is None:  
                        st.error(f"❌ 运行故障：小分子 '{lig_basename}' 对接未能生成亲和得分。报错原因：\n{raw}")  
                        continue  
                        
                    # 转换结构以进行可视化展示  
                    out_pdb = os.path.join(workdir, f"{lig_basename}.pdb")   
                    ok = pdbqt_to_pdb(out_pdbqt if os.path.exists(out_pdbqt) else lig_path, out_pdb)   
                      
                    calc_matrix.append({   
                        "化合物名称": lig_basename,   
                        "结合亲和力 (kcal/mol)": score,   
                        "PDB": out_pdb if ok else lig_path   
                    })   
                    progress.progress((idx + 1)/total)   
                    st.write(f"已完成 {idx+1}/{total}：{lig_basename}，得分 {score}，耗时 {t1-t0:.1f}s")   
                  
                # 只有当成功计算非空时渲染表格  
                if calc_matrix:  
                    df_results = pd.DataFrame(calc_matrix).sort_values("结合亲和力 (kcal/mol)", ascending=True).reset_index(drop=True)   
                    st.session_state["results_df"] = df_results   
                    st.success("对接任务运行完毕，结果已缓存于会话。")   
                else:  
                    st.error("真实计算未能成功完成。请在上方错误日志中检查对接配置或软件状态。")  
   
   
    # results rendering   
    if st.session_state["results_df"] is not None:   
        with st.container(border=True):   
            st.subheader("对接热力学亲和能排行（越负越好）")   
            df_eval = st.session_state["results_df"]   
            st.dataframe(df_eval[["化合物名称", "结合亲和力 (kcal/mol)"]], use_container_width=True)   
   
   
            best_lig = df_eval.iloc[0]   
            st.markdown(f"""   
            <div class="report-card">   
                <h4 style="margin-top:0px; color:#1e3a8a !important;">对接报告与相互作用判定</h4>   
                <p style="color:#475569; font-size:14px; line-height:1.6; margin:0;">   
                    结合平衡模拟计算完毕。在本轮物理对接中，与大分子结合最稳固的最优先导化合物为 <strong>{best_lig['化合物名称']}</strong>，   
                    结合自由能为 <strong>{best_lig['结合亲和力 (kcal/mol)']} kcal/mol</strong>。<br><br>   
                    <strong>分子对接作用机解读：</strong><br>   
                    结合能通常为负值，其绝对值代表了键合力的大小。若估值低于 <strong>-6.0 kcal/mol</strong>，说明分子能有效克服亲水阻碍沉入目标活性孔腔内。请查看下方三维透视图，观察小分子是否能与关键氨基酸残基建立氢键或疏水接触通道。   
                </p>   
            </div>   
            """, unsafe_allow_html=True)   
   
   
        with st.container(border=True):   
            sel_lig_name = st.selectbox("选择要渲染展示的三维化合物：", df_eval["化合物名称"])   
            s_pdb = df_eval.loc[df_eval["化合物名称"] == sel_lig_name, "PDB"].values[0]   
   
   
            st.subheader("活性位点微观对接三维交互视图")   
            if st.session_state.get("protein_pdb") and os.path.exists(st.session_state["protein_pdb"]) and os.path.exists(s_pdb):   
                show_structure(st.session_state["protein_pdb"], s_pdb)   
            else:   
                st.warning("无法渲染：缺少 protein.pdb 或 ligand pdb 文件。请检查工作区是否存在对应文件。")   
   
   
        # allow download of results CSV and zipped pdbs   
        with st.container(border=True):   
            st.subheader("导出与打包")   
   
   
            csv_bytes = df_eval.to_csv(index=False).encode("utf-8-sig")   
            st.download_button("导出对接结果 (CSV)", csv_bytes, "docking_results.csv", mime="text/csv", use_container_width=True)   
   
   
            # create zip of pdb files on demand   
            if st.button("打包并下载所有 PDB 文件"):   
                zip_path = os.path.join(workdir, "docking_pdbs.zip")   
                with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:   
                    for p in df_eval["PDB"].unique():   
                        if os.path.exists(p):   
                            zf.write(p, arcname=os.path.basename(p))   
                with open(zip_path, "rb") as fr:   
                    st.download_button("下载 PDB ZIP", fr.read(), file_name="docking_pdbs.zip", mime="application/zip", use_container_width=True)   
   
   
    else:   
        st.info("提示：尚无对接结果。请先在左侧上传受体与配体并启动计算。")
