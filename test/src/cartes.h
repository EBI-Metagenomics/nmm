#ifndef CARTES_H
#define CARTES_H

inline static int cartes_ipow(int base, int exp)
{
    int result = 1;
    for (;;) {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }

    return result;
}

inline static int cartes_nitems(int seq_len, int times)
{
    return cartes_ipow(seq_len, times);
}

inline static void cartes_item(char const* seq, int seq_len, int times, int idx, char* item)
{
    for (int i = 0; i < times; ++i) {
        item[i] = seq[(idx % cartes_ipow(seq_len, i + 1)) / cartes_ipow(seq_len, i)];
    }
}

#endif
