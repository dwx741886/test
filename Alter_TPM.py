# @Version   :1.0
# @Author    :杜伟轩
# @File      :Alter_TPM.py
# @Time      :2025/8/28 17:14
input_data = input("输入文件(csv):")
output_data = input("输出文件(csv):")
sample_content = int(input("样本量(整数):"))


def get_file():
    with open(input_data.strip(" "), "r") as f1:
        f1_list = f1.readlines()
        head = f1_list.pop(0).split(",")
        head.pop(-1)
        line_all = []
        for line in f1_list:
            l = line.replace("\n", "").split(",")
            for i in range(1, len(l)):
                l[i] = int(l[i])
            for j in range(1, sample_content + 1):
                a = l[j] / (l[sample_content + 1] / 1000)
                l.append(a)
            b = [l[0]]
            for k in range(sample_content + 2, 2 * sample_content + 2):
                b.append(l[k])
            line_all.append(b)
        return line_all, head


def get_sum(file1):
    sum_all = []
    for i in range(1, sample_content + 1):
        sample1 = []
        for line in file1:
            sample1.append(line[i])
        sum1 = 0
        for s in sample1:
            sum1 += s
        sum_all.append(sum1)
    return sum_all


def get_TPM(f1, l1, h1):
    ll = []
    for line in f1:
        l = [line[0]]
        for i in range(1, sample_content + 1):
            a = 1000000 * (line[i] / l1[i - 1])
            l.append(a)
        ll.append(l)
    with open(output_data.strip(" "), "w") as f:
        f.write(str(h1).strip("[").strip("]") + "\n")
        for line in ll:
            f.write(str(line).strip("[").strip("]") + "\n")


if __name__ == "__main__":
    f1, h1 = get_file()
    l1 = get_sum(f1)
    get_TPM(f1, l1, h1)
