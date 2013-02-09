public class MathUtils {
    public static long gcd(long a, long b) {
        if (b == 0) {
            return a;
        }
        return gcd(b, a % b);
    }

    public static long lcm(long a, long b) {
        return a * b / gcd(a, b);
    }

    public static long[] exgcd(long a, long b) {
        if (b == 0) {
            return new long[]{a,1,0};
        }
        long[] tmp = exgcd(b, a % b);
        return new long[]{tmp[0], tmp[2], (tmp[1] - (a/b) * tmp[2])};
    }

    public static int[] cycleFinding(int x0) {
        // f(x) is the function to be implemented
        int tortoise = f(x0), hare = f(f(x0));
        while (tortoise != hare) {
            tortoise = f(tortoise);
            hare = f(f(hare));
        }
        int[] ret = new int[]{0,0};
        hare = tortoise; tortoise = x0;
        while (tortoise != hare) {
            tortoise = f(tortoise);
            hare = f(hare);
            ret[0]++;
        }
        ret[1] = 1;
        hare = f(tortoise);
        while (tortoise != hare) {
            hare = f(hare);
            ret[1]++;
        }
        return ret;
    }

    /**
     * Euler's totient function (per query)
     **/
    public static long phiQuery(long N) {
        if(N == 1) {
            return 0;
        }
        long ret = N;
        for(long d = 2; d * d <= N; d++) if(N % d == 0) {
            ret -= ret / d;
            while(N % d == 0) {
                N /= d;
            }
        }
        if(N > 1) {
            ret -= ret / N;
        }
        return ret;
    }

    /**
     * Euler's totient function (pre-computation)
     **/
    final static int MAXN = (int)1e7;
    static int[] phi = new int[MAXN];
    public static void precomputePhi() {
        //IMPORTANT: Take care of the base cases and handle them alone wrt to the problem.
        phi[0] = 0; phi[1] = 0;
        for (int i = 2; i < MAXN; i++)  {
            phi[i] = i;
        }
        for (int i = 2; i < MAXN; i++) if (phi[i] == i) { // then i is a prime number
            for (int j = i; j < MAXN; j += i) {
                phi[j] = phi[j] / i * (i - 1);  // divide by i then multiply by i-1
            }
        }
    }

    /**
     * Pre-computing of Prime factors
     **/
    //final static int MAXN = (int)1e7; //declared before precomputePhi() function
    static int primeFactor[] = new int[MAXN]; //holds a single prime factor for each number,initialized with 0
    private static void calculatePrimeFactors() {
        for (int i = 2; i * i < MAXN; i++) {
            // we didn't iterate over the whole interval cause we only need to
            // know a single prime factor for each composite
            if (primeFactor[i] == 0) {
                for (int j = i * i; j < MAXN; j += i)
                    primeFactor[j] = i;
            }
        }
        for (int i = 2; i < MAXN; i++) if (primeFactor[i] == 0) {
                primeFactor[i] = i; // sets every prime factor of a prime number as the prime it self
        }
        int[] nPrimeFactors = new int[MAXN]; // we'll have the result here
        nPrimeFactors[0] = nPrimeFactors[1] = 0;
        for (int i = 2; i < MAXN; i++) {
            nPrimeFactors[i] = 1 + nPrimeFactors[i / primeFactor[i]]; // Like a DP ;)
        }
    }

    /**
     * Gauss Jordan method for solving a system of linear equations
     **/
    double[][] A;    // take care array A is updated after calling gauss()
    double[] gauss() {
        double[] sol = new double[A.length];
        int[] order = new int[A.length];
        Arrays.fill(order, -1);	
        for (int column = 0; column < A[0].length-1; column++) {
            //get max value in this column and not visited
            int index = -1;
            double value = - Double.MAX_VALUE;
            for (int row = 0; row < A.length; row++) {
                if (order[row]==-1 && A[row][column] > value &&
                        !equal(A[row][column],0.0)) {
                    value = A[row][column];
                    index = row;
                }
            }
            order[index] = column;
            for (int i = 0; i < A.length; i++) {
                if (i == index) continue;
                for (int j = column+1; j < A[i].length; j++) {
                    double val = A[i][j] * A[index][column];
                    val -= A[i][column] * A[index][j];
                    val /= A[index][column];
                    A[i][j] = val;
                }
            }		
            for (int i = column + 1; i < A[index].length; i++) {
                A[index][i] /= -A[index][column];
            }
        }
        for (int i = 0; i < sol.length; i++) {
            sol[order[i]] = A[i][A[i].length - 1];
        }
        return sol;
    }
}
