# 🧬 生物序列翻译工具

一个简单易用的生物序列翻译工具，支持从NCBI下载DNA/RNA序列并翻译为蛋白质序列。

## ✨ 功能特性

- ✅ **从NCBI下载序列**: 支持输入Accession ID（如：NM_001301717）自动下载序列
- ✅ **序列清理**: 自动去除空格和非法字符
- ✅ **蛋白质翻译**: 将DNA/RNA序列翻译为蛋白质序列
- ✅ **阅读框支持**: 支持设置起始位点和阅读框（1、2、3）
- ✅ **密码子对照表**: 生成详细的密码子→氨基酸对照表
- ✅ **结果导出**: 支持下载CSV格式的密码子表
- ✅ **双版本支持**: 提供网页版（Streamlit）和命令行版本

## 🚀 快速开始

### 环境要求

- Python 3.7+
- 网络连接（用于NCBI数据下载）

### 安装依赖

#### 方法1：使用虚拟环境（推荐）

```bash
# 创建虚拟环境
python -m venv venv

# 激活虚拟环境
# Windows:
.\venv\Scripts\activate
# Linux/Mac:
source venv/bin/activate

# 安装依赖
pip install -r requirements.txt
```

#### 方法2：全局安装

```bash
pip install -r requirements.txt
```

### 使用方式

#### 1. 网页版（推荐）

启动Streamlit应用：

```bash
# 使用虚拟环境
.\venv\Scripts\activate
streamlit run app.py

# 或使用批处理文件（Windows）
run_app_venv.bat
```

然后在浏览器中打开显示的地址（通常是 http://localhost:8501）

#### 2. 命令行版本

```bash
# 使用NCBI Accession ID
python cli_translator.py -a NM_001301717 -e your-email@example.com

# 直接输入序列
python cli_translator.py -s "ATGGCCGAA..." -e your-email@example.com

# 设置起始位置和阅读框
python cli_translator.py -a NM_001301717 -e your-email@example.com -p 10 -f 2

# 自定义输出文件名
python cli_translator.py -a NM_001301717 -e your-email@example.com -o my_result
```

## 📖 使用说明

### 网页版功能

1. **序列输入**
   - 选择输入方式：NCBI Accession ID 或 直接输入序列
   - 输入Accession ID后点击"下载序列"
   - 或直接在文本框中粘贴DNA/RNA序列

2. **配置选项**（侧边栏）
   - 选择阅读框（1、2、3）
   - 设置起始位置（从0开始）

3. **结果查看**
   - 序列信息：长度、GC含量、序列预览
   - 翻译结果：蛋白质序列、氨基酸组成
   - 密码子对照表：详细的位置、密码子、氨基酸对应关系

4. **结果导出**
   - 点击"下载密码子表"按钮保存CSV文件

### 命令行参数

| 参数 | 简写 | 说明 | 必需 |
|------|------|------|------|
| `--accession` | `-a` | NCBI Accession ID | 与 `-s` 二选一 |
| `--sequence` | `-s` | 直接输入DNA/RNA序列 | 与 `-a` 二选一 |
| `--email` | `-e` | NCBI邮箱地址 | ✅ |
| `--start` | `-p` | 起始位置（默认：0） | ❌ |
| `--frame` | `-f` | 阅读框（1/2/3，默认：1） | ❌ |
| `--output` | `-o` | 输出文件前缀（默认：result） | ❌ |

## 📁 输出文件

### 网页版
- 密码子对照表：CSV格式，包含位置、密码子、氨基酸信息

### 命令行版
- `{prefix}_protein.fasta`: 翻译后的蛋白质序列（FASTA格式）
- `{prefix}_codon_table.csv`: 密码子对照表（CSV格式）

## 🔧 技术实现

- **后端**: Python 3.7+
- **生物信息学**: Biopython
- **网页界面**: Streamlit
- **数据处理**: Pandas
- **网络请求**: Requests

## 📝 注意事项

1. **NCBI邮箱**: 使用NCBI服务需要提供有效的邮箱地址
2. **网络连接**: 下载序列需要稳定的网络连接
3. **序列格式**: 支持DNA（ATCG）和RNA（AUCG）序列
4. **阅读框**: 阅读框1从第1位开始，阅读框2从第2位开始，阅读框3从第3位开始

## 🤝 贡献

欢迎提交Issue和Pull Request来改进这个工具！

## 📄 许可证

MIT License

## 🔗 相关链接

- [Biopython官方文档](https://biopython.org/)
- [NCBI数据库](https://www.ncbi.nlm.nih.gov/)
- [Streamlit文档](https://docs.streamlit.io/)