# @Version   :1.0
# @Author    :杜伟轩
# @File      :merge_left.py
# @Time      :2025/8/16 22:42
# import pandas as pd  # 给一个模块别名，方便后续使用

def align_file(key, file):
    with open(key, "r") as f1:
        with open(file, "r") as f2:
            f1_list = f1.readlines()
            f2_list = f2.readlines()
            f1_list.pop(0)
            f2_list.pop(0)
            key_list = []
            for i in f1_list:
                key_list.append(i.replace("\n", ""))
            obj_list = []
            for j in f2_list:
                data = j.replace("\n", "").split(",")
                obj_list.append(data)
            re_list = []
            for k in key_list:
                for l in obj_list:
                    if l[0] == k:
                        re_list.append(l)
            return re_list

def write_file(file, lists):
    with open(file, "w") as f:
        for line in lists:
            f.write(str(line).strip("[").strip("]")+"\n")

if __name__== "__main__":
    input_file1 = input("请输入文件[key]（csv）:")
    input_file2 = input("请输入文件（列名要有[key]，（csv））:")
    output = input("输出文件（csv）:")
    aligned_file = align_file(input_file1,input_file2)
    write_file(output, aligned_file)



