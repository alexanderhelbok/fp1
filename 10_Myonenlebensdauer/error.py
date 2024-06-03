import uncertainties as unc
import uncertainties.unumpy as unp

def contributions(var, rel=True, precision=2):
    if rel:
        for (name, error) in var.error_components().items():
            print("{}: {} %".format(name.tag, round(error ** 2 / var.s ** 2 * 100, precision)))
    else:
        for (name, error) in var.error_components().items():
            print("{}: {}".format(name.tag, round(error, precision)))

texp = unc.ufloat(3.3, 0.1, "time")

b = unc.ufloat(2.00, 0.03, "tau")
c = unc.ufloat(17.9, 0.5, "offset")
a = unc.ufloat(24100, 200, "amplitude")

def model1(x, a, b):
    return a/(10*b) * unp.exp(-x / b)


def model2(x, a, b, c):
    return a/(10*b) * unp.exp(-x / b) + c

out = model1(texp, a, b)/model2(texp, a, b, c)

print(f"{out:.1uS}")
contributions(out, precision=1)