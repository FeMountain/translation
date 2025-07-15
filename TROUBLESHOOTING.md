# 🔧 故障排除指南

## Linter 导入错误解决方案

如果您在VS Code或其他IDE中看到类似以下的导入错误：
```
Import "streamlit" could not be resolved
Import "pandas" could not be resolved  
Import "Bio" could not be resolved
```

### 解决方案

#### 方法1：使用虚拟环境（推荐）

1. **激活虚拟环境**：
   ```bash
   # Windows
   .\venv\Scripts\activate
   
   # Linux/Mac
   source venv/bin/activate
   ```

2. **在VS Code中选择正确的Python解释器**：
   - 按 `Ctrl+Shift+P`
   - 输入 "Python: Select Interpreter"
   - 选择 `./venv/Scripts/python.exe`

3. **重启VS Code**：
   - 关闭VS Code
   - 重新打开项目
   - Linter错误应该消失

#### 方法2：手动指定Python路径

在VS Code的settings.json中添加：
```json
{
    "python.defaultInterpreterPath": "./venv/Scripts/python.exe",
    "python.analysis.extraPaths": ["./venv/Lib/site-packages"]
}
```

#### 方法3：全局安装依赖

如果不想使用虚拟环境：
```bash
pip install -r requirements.txt
```

### 验证安装

运行以下命令验证所有包是否正确安装：
```bash
python -c "import streamlit, pandas, Bio; print('✅ 所有包导入成功')"
```

### 启动应用

#### 网页版
```bash
# 使用虚拟环境
.\venv\Scripts\activate
streamlit run app.py

# 或使用批处理文件
run_app_venv.bat
```

#### 命令行版
```bash
# 使用虚拟环境
.\venv\Scripts\activate
python cli_translator.py -s "ATGGCCGAA" -e your-email@example.com
```

## 常见问题

### Q: 为什么会有linter错误？
A: Linter错误通常是因为IDE无法找到已安装的包。使用虚拟环境可以解决这个问题。

### Q: 应用能正常运行但IDE显示错误？
A: 这是正常的，说明包已正确安装，只是IDE配置问题。按照上述方法配置IDE即可。

### Q: 如何确认包已正确安装？
A: 运行测试脚本：`python test_example.py` 