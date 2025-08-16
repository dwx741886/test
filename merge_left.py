# @Version   :1.0
# @Author    :杜伟轩
# @File      :merge_left.py
# @Time      :2025/8/16 22:42
import pandas as pd  # 给一个模块别名，方便后续使用

input_file1 = input("请输入文件[key]（csv）:")
input_file2 = input("请输入文件（列名要有[key]，（csv））:")
output = input("输出文件（csv）:")
# pd.set_option('display.max_rows', None)  # 显示所有列，可忽略，不影响
# pd.set_option('display.max_columns', None)  # 显示所有行，可忽略，不影响
df1 = pd.read_csv(input_file1, names=["key"], header=0)  # id第一列空白，去掉列名，用names定义基因列名字
# df1.head()  # 查看df1数据
df2 = pd.read_csv(input_file2,header=0)  # 需有列名，用header定义，注意列名一定有一个是“key”，需检查df1和df2存在相同列名
# df2.head()
df_merge = pd.merge(left=df1, right=df2, on='key', how='left')
# print("/gss1/home/chenjf01/", df_merge)  # 只是打印展示
df_merge.to_csv("/Users/ackee_z/zxk/study/R/chc_count_TY.csv")
