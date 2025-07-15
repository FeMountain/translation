@echo off
echo 🧬 启动生物序列翻译工具 (虚拟环境版本)
echo.
echo 激活虚拟环境...
call venv\Scripts\activate.bat
echo.
echo 检查Python环境...
python --version
echo.
echo 检查依赖包...
python -c "import streamlit, pandas, Bio; print('✅ 所有依赖包已安装')"
echo.
echo 启动Streamlit应用...
streamlit run app.py
pause 