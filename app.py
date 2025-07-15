import streamlit as st
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import re
import io

# 设置页面配置
st.set_page_config(
    page_title="生物序列翻译工具",
    page_icon="🧬",
    layout="wide"
)

# 设置NCBI邮箱（必需）
Entrez.email = "your-email@example.com"

def clean_sequence(sequence):
    """清理序列，去除空格和非法字符"""
    # 去除空格、换行符等
    cleaned = re.sub(r'\s+', '', sequence)
    # 只保留DNA/RNA字符
    cleaned = re.sub(r'[^ATCGUatcgu]', '', cleaned)
    return cleaned.upper()

def download_sequence(accession_id):
    """从NCBI下载序列"""
    try:
        # 搜索序列
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq), record.description
    except Exception as e:
        st.error(f"下载序列失败: {str(e)}")
        return None, None

def translate_sequence(sequence, start_position=0, reading_frame=1):
    """翻译DNA/RNA序列为蛋白质"""
    try:
        # 创建Seq对象
        seq_obj = Seq(sequence)
        
        # 根据阅读框调整起始位置
        if reading_frame == 2:
            start_position += 1
        elif reading_frame == 3:
            start_position += 2
        
        # 从指定位置开始翻译
        if start_position > 0:
            seq_obj = seq_obj[start_position:]
        
        # 翻译序列
        protein_seq = seq_obj.translate(to_stop=False)
        
        return str(protein_seq)
    except Exception as e:
        st.error(f"翻译失败: {str(e)}")
        return None

def get_codon_table(sequence, start_position=0, reading_frame=1):
    """生成密码子到氨基酸的对照表"""
    try:
        # 调整起始位置
        if reading_frame == 2:
            start_position += 1
        elif reading_frame == 3:
            start_position += 2
        
        # 从指定位置开始
        if start_position > 0:
            sequence = sequence[start_position:]
        
        # 确保序列长度是3的倍数
        if len(sequence) % 3 != 0:
            sequence = sequence[:-(len(sequence) % 3)]
        
        codons = []
        amino_acids = []
        positions = []
        
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3]
            if len(codon) == 3:
                # 创建临时Seq对象进行翻译
                temp_seq = Seq(codon)
                aa = str(temp_seq.translate(to_stop=False))
                
                codons.append(codon)
                amino_acids.append(aa)
                positions.append(f"{start_position + i + 1}-{start_position + i + 3}")
        
        return pd.DataFrame({
            '位置': positions,
            '密码子': codons,
            '氨基酸': amino_acids
        })
    except Exception as e:
        st.error(f"生成密码子表失败: {str(e)}")
        return None

def main():
    st.title("🧬 生物序列翻译工具")
    st.markdown("---")
    
    # 侧边栏配置
    st.sidebar.header("⚙️ 配置选项")
    
    # 输入方式选择
    input_method = st.sidebar.radio(
        "选择输入方式",
        ["NCBI Accession ID", "直接输入序列"]
    )
    
    # 阅读框设置
    reading_frame = st.sidebar.selectbox(
        "选择阅读框",
        [1, 2, 3],
        format_func=lambda x: f"阅读框 {x}"
    )
    
    # 起始位置设置
    start_position = st.sidebar.number_input(
        "起始位置 (从0开始)",
        min_value=0,
        value=0,
        step=1
    )
    
    # 主界面
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.subheader("📥 序列输入")
        
        if input_method == "NCBI Accession ID":
            accession_id = st.text_input(
                "输入NCBI Accession ID",
                placeholder="例如: NM_001301717"
            )
            
            if st.button("🔍 下载序列", type="primary"):
                if accession_id:
                    with st.spinner("正在从NCBI下载序列..."):
                        sequence, description = download_sequence(accession_id)
                        if sequence:
                            st.success(f"✅ 成功下载序列: {description}")
                            st.session_state.sequence = sequence
                            st.session_state.description = description
                else:
                    st.warning("请输入Accession ID")
        else:
            sequence_input = st.text_area(
                "直接输入DNA/RNA序列",
                placeholder="输入ATCG序列...",
                height=200
            )
            
            if st.button("🧹 清理并处理序列", type="primary"):
                if sequence_input:
                    cleaned_sequence = clean_sequence(sequence_input)
                    if cleaned_sequence:
                        st.success(f"✅ 序列清理完成，长度: {len(cleaned_sequence)} bp")
                        st.session_state.sequence = cleaned_sequence
                        st.session_state.description = "用户输入序列"
                else:
                    st.warning("请输入序列")
    
    with col2:
        st.subheader("📊 序列信息")
        if 'sequence' in st.session_state:
            sequence = st.session_state.sequence
            st.metric("序列长度", f"{len(sequence)} bp")
            st.metric("GC含量", f"{((sequence.count('G') + sequence.count('C')) / len(sequence) * 100):.1f}%")
            
            # 显示序列预览
            st.text("序列预览:")
            st.code(sequence[:100] + "..." if len(sequence) > 100 else sequence)
    
    # 翻译结果
    if 'sequence' in st.session_state:
        st.markdown("---")
        st.subheader("🔄 翻译结果")
        
        # 翻译序列
        protein_sequence = translate_sequence(
            st.session_state.sequence,
            start_position,
            reading_frame
        )
        
        if protein_sequence:
            col1, col2 = st.columns(2)
            
            with col1:
                st.metric("蛋白质长度", f"{len(protein_sequence)} aa")
                st.text("蛋白质序列:")
                st.code(protein_sequence)
            
            with col2:
                # 统计氨基酸组成
                aa_counts = {}
                for aa in protein_sequence:
                    if aa != '*':  # 排除终止密码子
                        aa_counts[aa] = aa_counts.get(aa, 0) + 1
                
                if aa_counts:
                    st.text("氨基酸组成:")
                    for aa, count in sorted(aa_counts.items()):
                        st.text(f"{aa}: {count}")
            
            # 密码子对照表
            st.markdown("---")
            st.subheader("📋 密码子对照表")
            
            codon_table = get_codon_table(
                st.session_state.sequence,
                start_position,
                reading_frame
            )
            
            if codon_table is not None:
                st.dataframe(codon_table, use_container_width=True)
                
                # 下载选项
                csv = codon_table.to_csv(index=False, encoding='utf-8-sig')
                st.download_button(
                    label="📥 下载密码子表 (CSV)",
                    data=csv,
                    file_name=f"codon_table_{start_position}_{reading_frame}.csv",
                    mime="text/csv"
                )

if __name__ == "__main__":
    main() 