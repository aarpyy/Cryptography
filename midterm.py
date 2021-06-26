from crypto_functions import *

power = 2
p = RandomPrime(10000)

lst = [x for x in range(p)]
nums = {}
for e in lst:
    nums[e] = 0
for e in range(1, p):
    k = pow(e, power, p)
    if k in nums:
        nums[k] += 1
    else:
        nums[k] = 1

nums2 = {}
for n in nums:
    r = nums[n]
    if r in nums2:
        nums2[r] += 1
    else:
        nums2[r] = 1


#print(nums2)
print("You see:", end=" ")
for n in nums2:
    k = nums2[n]
    print(f"{k} values {n} times;", end=" ")
print()
print(f"GCD: {GCD(power, p-1)}, prime: {p}")
# power = 9
# p = 6637
#
# nums = {}
# for a in range(1, p):
#     k = pow(a, power, p)
#     if k in nums:
#         nums[k].append(a)
#     else:
#         nums[k] = [a]
# print(nums)