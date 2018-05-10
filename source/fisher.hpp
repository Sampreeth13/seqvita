#include <boost/math/distributions/hypergeometric.hpp>
#include <algorithm> // for min and max

using namespace boost::math;
using namespace std;

double fisher_test(unsigned a, unsigned b, unsigned c, unsigned d) {
		unsigned N = a + b + c + d;
		unsigned r = a + c;
		unsigned n = c + d;
		unsigned max_for_k = min(r, n);
		unsigned min_for_k = (unsigned)max(0, int(r + n - N));
		hypergeometric_distribution<> hgd(r, n, N);
		double tmp_p = pdf(hgd,c);
		if((a*d) >= (b*c)){
			unsigned int low = min(b,c);
			for(int i=0;i<low;i++){
				a++;b--;c--;d++;
				tmp_p += pdf(hgd,c);
			}
		}
		else{
			unsigned int low = min(a,d);
			for(int i=0;i<low;i++){
				a--;b++;c++;d--;
				tmp_p += pdf(hgd,c);
			}
		}
		return tmp_p;
}
