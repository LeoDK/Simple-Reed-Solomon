class RS:

    def __init__ (self, d, k, t):
        self.F = GF(2^d, 'a')
        self.X = PolynomialRing(self.F, 'X').gen()
        self.alpha = self.F.gen()
        self.k = k
        self.t = t
        self.n = k+2*t
        self.d = d
        self.G = self.X^0
        for i in range (1, 2*self.t+1):
            self.G *= (self.X - self.alpha^i)

    def encode (self, m):
        # [m] is a list of elements between 0 and 2^self.d - 1
        P = 0 * self.X
        for i in range(len(m)):
            P += self.F.fetch_int(m[i]) * self.X^i
        return P*self.G

    @staticmethod
    def int_repr (P):
        # Return P in integer representation
        return [ coef.integer_representation() for coef in P.list() ]

    def _find_power (self, elem):
        # Find i such that elem = self.alpha^i
        curr = self.alpha^0
        for i in range(self.F.order()+1):
            if curr == elem:
                return i
            curr *= self.alpha

    def correct (self, D):
        # [D] is a polynomial in self.F[X] of degree less than k+2t

        # First step : finding the number of errors and find Lambda
        S = vector(self.F, [D(self.alpha^i) for i in range(2*self.t+1)])

        if S == [0]*(2*self.t+1): # The received message is perfect
            return (0,D)

        M = matrix(self.F, self.t)
        for i in range(self.t):
            for j in range(self.t):
                M[i,j] = S[self.t + i - j]

        # Reducing until we find a non zero discriminant
        nu = self.t
        while (M.determinant() == 0):
            M = M[:-1,1:]
            nu -= 1
        lambdas = M^-1 * S[nu+1:2*nu+1]

        # Second step : find the x_r
        Lambda = self.X^0
        for i in range(len(lambdas)):
            Lambda += lambdas[i] * self.X^(i+1)

        roots = Lambda.roots()
        x = [root[0]^-1 for root in roots]

        # Third step : find the places where are located the errors
        location = [self._find_power(x_r) for x_r in x]

        # Fourth step : find the y_r (errors at each place)
        N = matrix(self.F, nu)
        for i in range(nu):
            for j in range(nu):
                N[i,j] = x[j]^(i+1)

        y = N^-1 * S[1:nu+1]

        # Fifth step : recover the real message
        E = 0 * self.X
        for (i_r,y_r) in zip(location, y):
            E += y_r * self.X^i_r
        C = D+E

        return (nu, C)

    def decode (self, D):
        nu, C = self.correct(D)
        P = C//self.G
        m = [coef.integer_representation() for coef in P.list()]
        return (nu, m)
