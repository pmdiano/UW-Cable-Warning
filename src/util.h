#ifndef __UTIL_H__
#define __UTIL_H__

// Get the sum of a 2-D array
template<typename T>
T sumImage(T **BWBscope, int M, int N)
{
    T s = static_cast<T>(0);

    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
        {
            s += BWBscope[i][j];
        }
    }

    return s;
}

int clip(double value, int low, int high)
{
    int ret = static_cast<int>(value + 0.5);

    if (ret < low)
    {
        return low;
    }
    if (ret > high)
    {
        return high;
    }
    return ret;
}

#endif
