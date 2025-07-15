#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
测试示例：验证翻译功能
"""

from Bio.Seq import Seq

def test_translation():
    """测试基本的翻译功能"""
    print("🧪 测试翻译功能")
    
    # 测试序列（一个简单的DNA序列）
    test_sequence = "ATGGCCGAA"
    
    print(f"测试序列: {test_sequence}")
    
    # 创建Seq对象并翻译
    seq_obj = Seq(test_sequence)
    protein = seq_obj.translate()
    
    print(f"翻译结果: {protein}")
    print(f"预期结果: MAE")
    
    if str(protein) == "MAE":
        print("✅ 翻译功能正常")
    else:
        print("❌ 翻译功能异常")

if __name__ == "__main__":
    test_translation() 