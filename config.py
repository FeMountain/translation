# -*- coding: utf-8 -*-
"""
配置文件
"""

# NCBI设置
NCBI_EMAIL = "your-email@example.com"  # 请替换为您的邮箱

# 应用设置
APP_TITLE = "🧬 生物序列翻译工具"
APP_ICON = "🧬"

# 默认参数
DEFAULT_READING_FRAME = 1
DEFAULT_START_POSITION = 0
DEFAULT_OUTPUT_PREFIX = "result"

# 序列设置
MAX_SEQUENCE_LENGTH = 100000  # 最大序列长度
MIN_SEQUENCE_LENGTH = 3       # 最小序列长度

# 文件设置
SUPPORTED_FORMATS = ['.fasta', '.fa', '.txt', '.csv'] 