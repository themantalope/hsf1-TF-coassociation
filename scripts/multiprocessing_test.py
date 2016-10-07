import multiprocessing


def untuplefoo(tup):
    return foo(tup[0], tup[1])

def foo(bar, baz):
    return bar + baz


inputs = [(i, i+1) for i in range(0, 20,2)]


# print inputs

p = multiprocessing.Pool(multiprocessing.cpu_count())

print p.map(untuplefoo, inputs)
a = p.map(untuplefoo, inputs)
print a
# for i in p.imap(untuplefoo, inputs):
#     print i