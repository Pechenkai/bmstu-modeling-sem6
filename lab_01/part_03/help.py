import math

def euler_one_step(x, U, step):
    f_val = step * x + U * step * U * U
    big_result = U + f_val
    return big_result


def find_const(x, U):
    return 1 / (2 * U * U) + x

def simpson(x_left, x_right, y_left, y_right):
    dx = x_right - x_left
    x_mid = (x_left + x_right) / 2.0
    y_mid = (y_left + y_right) / 2.0

    f_left = x_left + y_left ** 3
    f_right = x_right + y_right ** 3
    f_mid = x_mid + y_mid ** 3

    return (dx / 6.0) * (f_left + 4.0 * f_mid + f_right)


def picard(x0, h, n, num_iter):
    res = []

    prev = [0.0] * n
    x = x0
    x1 = x0 + h

    for i in range(1, num_iter):
        curr = [0.0] * n

        x_left = x
        x_right = x1

        y_left = prev[i - 1]
        y_right = prev[i]

        integral = simpson(x_left, x_right, y_left, y_right)

        curr[i] = curr[i - 1] + integral

        res.append(curr)
        prev = curr
        x_left += h
        x_right += h

    return res


def adaptive_euler(eps=1e-4, h_init=1, x_max_limit=1.65):
    x = 0.0
    U = 0.0
    h = h_init
    h_limit = 1e-323

    f = 1
    f10 = 0
    U_prev = 0.0

    while True:
        U_big = euler_one_step(x, U, h)

        half = h / 2.0
        U_half_1 = euler_one_step(x, U, half)
        U_small = euler_one_step(x + half, U_half_1, half)

        denominator = max(1.0, abs(U_small))
        local_error = abs(U_big - U_small) / denominator

        if local_error > eps:
            h /= 2.0
            if h < h_limit:
                print("Шаг стал слишком малым")
                return x, U_small, h
            continue
        else:
            x += h
            U = U_small

            if f10 < 10 and U > 10 ** f10:
                print(f" x={x:.20E}\n c={find_const(x, U):.80E}")
                f10 += 1

            if U > 10 ** (f * 10):
                print(f"x={x:.20E}, U={U:.8E}, h={h:.8E}")
                print(f"c={find_const(x, U):.80E}")
                f += 1

            if math.isinf(U):
                return x - h, U_prev, h * 2

            if x > x_max_limit:
                return x_max_limit, U, h

            U_prev = U

# def task3():
#     x0 = 0
#     x_max = 1.707246    # Рунге Кутта 4 порядка
#     # x_max = 1.823675    # Эйлер
#     N = 100
#     h = (x_max - x0) / N
#     Xvals = []
#     Yvals_first = []
#     Yvals_second = []
#     Yvals_third = []
#     Yvals_forth = []
#     Yvals_euler = euler_base(0, 0, h, N, func_3)
#     Yvals_runge_kutta = runge_kutta4(0, 0, h, N, func_3)
#     for i in range(N):
#         x = x0 + i * h
#         Xvals.append(x)
#         Yvals_first.append(pikar_first_3(x))
#         Yvals_second.append(pikar_second_3(x))
#         Yvals_third.append(pikar_third_3(x))
#         Yvals_forth.append(pikar_forth_3(x))
#     table_print(Xvals, Yvals_first, Yvals_second, Yvals_third, Yvals_forth, Yvals_euler, Yvals_runge_kutta)
#     plt.plot(Xvals, Yvals_euler, label='Эйлер')
#     plt.plot(Xvals, Yvals_runge_kutta, label='Рунге-Кутта 4')
#     plt.plot(Xvals, Yvals_first, label='1 приближение')
#     plt.plot(Xvals, Yvals_second, label='2 приближение')
#     plt.plot(Xvals, Yvals_third, label='3 приближение')
#     plt.plot(Xvals, Yvals_forth, label='4 приближение')
#     plt.legend()
#     plt.show()
#
# def table_print(x, y1, y2, y3, y4, y5, y6):
#     print("---------------------------------------------------------------------------------------------------------")
#     print("|      x      |    Пикар 1   |    Пикар 2   |    Пикар 3   |    Пикар 4   |     Эйлер    | Рунге-Кутта4 |")
#     print("---------------------------------------------------------------------------------------------------------")
#     for i in range(len(x)):
#         print("|{:^14.5f}|{:^14.5f}|{:^14.5f}|{:^14.5f}|{:^14.5f}|{:^14.5f}|{:^14.5f}|".format(
#             x[i], y1[i], y2[i], y3[i], y4[i], y5[i], y6[i]))
#     print("---------------------------------------------------------------------------------------------------------")
#
#
# def pikar_first_3(x):
#     return x ** 2 / 2
#
# def pikar_second_3(x):
#     return x ** 7 / 56 + x ** 2 / 2
#
# def pikar_third_3(x):
#     return (x ** 22 / 3863552
#             + 3 * x ** 17 / 106624
#             + x ** 12 / 896
#             + x ** 7 / 56
#             + x ** 2 / 2)
#
# def pikar_forth_3(x):
#     return (x ** 67 / 3863981943017739124736
#             + 9 * x ** 62 / 98677964914244452352
#             + 1081 * x ** 57 / 73440052228804050944
#             + 76623 * x ** 52 / 53388985337387155456
#             + 310503 * x ** 47 / 3244062457475366912
#             + (1069 * x ** 42) / 224099368435712
#             + (790667 * x ** 37) / 4145838316060672
#             + (11871 * x ** 32) / 1883187970048
#             + (361 * x ** 27) / 2101772288
#             + (47 * x ** 22) / 11941888
#             + (33 * x ** 17) / 426496
#             + (x ** 12) / 896
#             + (x ** 7) / 56
#             + (x ** 2) / 2)
#
# pub fn find_max_x_with_tolerance<F>(
#     f: F,
#     x0: f64,
#     u0: f64,
#     h: f64,
#     eps: f64,
# ) -> f64
# where
#     F: Fn(f64, f64) -> f64,
# {
#     let mut x_current = x0;
#     let mut u_big = u0;
#     let mut u_small = u0;
#
#     loop {
#         let f_val = f(x_current, u_big);
#         let u_big_next = u_big + h * f_val;
#         let x_next = x_current + h;
#
#         let half = h / 2.0;
#         let u_temp = u_small + half * f(x_current, u_small);
#         let x_mid = x_current + half;
#         let u_temp = u_temp + half * f(x_mid, u_temp);
#
#         let u_small_next = u_temp;
#
#         let denom = 1.0_f64.max(u_small_next.abs());
#         let rel_err = (u_big_next - u_small_next).abs() / denom;
#
#         if rel_err > eps {
#             return x_current;
#         }
#
#         x_current = x_next;
#         u_big = u_big_next;
#         u_small = u_small_next;
#     }
# }