#ifndef IDDH
#define IDDH

#include<vector>
#include<map>
#include<list>
#include<stdlib.h>
#include<math.h>
#include<gmp.h>
#include<iostream>

#ifndef TESTING
#define WORD_OPT
#endif

#define IDD_WORD_SIZE 5   // sizeof(WORD) = 2^5 = 32

// Status
#define idd_fit_ulong 1


using namespace std;
class IDDContent;
class uIDD;
class sIDD;


class IDDC  // Protected container
{
  private:
  IDDContent *a;
  
  public:
  IDDC();
  IDDC(const IDDC &s);
  IDDC(IDDContent* _a);
  ~IDDC();
  IDDC& operator = (const IDDC &m);
  IDDC& operator = (IDDContent *m);
  inline bool operator == (const IDDC &m) const {return a==m.a;}
  inline bool operator == (IDDContent *b) const {return a==b;}    
  friend bool operator == (const IDDContent *n, const IDDC &m);
  inline bool operator != (const IDDC &m) const {return a!=m.a;}  
  friend bool operator != (const IDDContent *n, const IDDC &m);
  inline IDDContent* operator ->() const {return a;}
  inline IDDContent* operator *() const {return a;}  
  inline operator unsigned int () const {return (unsigned int)a;}
  
  friend void idd_print(const uIDD &n);
  friend void idd_print(const sIDD &n);  
};
inline bool operator == (const IDDContent *n, const IDDC &m) {return n==m.a;}
inline bool operator != (const IDDContent *n, const IDDC &m) {return n!=m.a;}


//----------------------------------------------------------------------------//
//                                                                            //
//                               IDDC Container                                //
//                                                                            //
//----------------------------------------------------------------------------//
class IDDContent   
{
  private:
  // Hashtable
  template <class key, typename data = IDDContent*>
  class hash
  {
    private:
    vector<map<key, data>* > T;
    unsigned int n_items;
    unsigned int init_size;
    unsigned int max_size;
    bool (*onResizing)();
    
    int next_prime(int n);
    
    public:  
    hash(unsigned int size, bool (*resizing)() = NULL, unsigned int max = 0); 
    ~hash();   
    bool exists(const key &k);
    data find(const key &k);
    void insert(const key &k, const data &d);
    void resize(unsigned int size);
    void remove(const key &k);
    bool clean();
    void clear();
    inline unsigned int items_count() {return n_items;}
    inline unsigned int size() {return T.size();}
    inline unsigned int memory();
    void show_stats();
    friend class IDDContent;
  };
  class IDDC1
  {
    public:
    const IDDContent* n;
    IDDC1(const IDDContent* _n):n(_n) {}
    inline unsigned int hash() const {return (unsigned int)(n)^((unsigned int)(n)<<16);}
    inline bool operator < (const IDDC1 &m) const {return n<m.n;}
  };
  class IDDC2
  {
    public:
    const IDDContent *a, *b;
    IDDC2(const IDDContent* _a, const IDDContent* _b):a(_a),b(_b) {}
    inline unsigned int hash() const {return ((unsigned int)(a)^(unsigned int)(b))+((unsigned int)(a)<<16)-((unsigned int)(b)>>16);}
    inline bool operator < (const IDDC2 &m) const {return a<m.a || (a==m.a && b<m.b);}    
  };
  class IDDC3
  {
    public:
    const IDDContent *g, *p, *d;
    IDDC3(const IDDContent* _g, const IDDContent* _p, const IDDContent* _d):g(_g),p(_p),d(_d) {}
    inline unsigned int hash() const {return ((unsigned int)(g)^(unsigned int)(d))-(unsigned int)(p)+((unsigned int)(g)>>16)-((unsigned int)(d)<<16);}
    inline bool operator < (const IDDC3 &m) const {return g<m.g || (g==m.g && p<m.p) || (g==m.g && p==m.p && d<m.d);}    
  };  

  
  // Trichotomy polongers
  IDDContent* g;
  IDDContent* p;
  IDDContent* d;
  
  inline bool ul_fit() const {return p==NULL;}
  inline unsigned int val() const {return (unsigned int)g;}
  unsigned int get_p() const
  {
    unsigned int q;
    for(q = IDD_WORD_SIZE-1;
       (val()>>(1<<q))==0 && q!=0;
        q--) {}
    return q;
  }
  unsigned int get_g(const unsigned int &lp) const {return val()&((unsigned int)(1<<(1<<lp))-1);} 
  unsigned int get_d(const unsigned int &lp) const {return val()>>(1<<lp);}   
  
  // Is root?
  unsigned int n_shared;
  static bool clean_Ht();

  // Global memo hash-tables
  static hash<IDDC3>  Ht;
  static hash<IDDC1>  Hx;
  static IDDC zero;
  static IDDC one;  
  static IDDC Xws;
    
  // internal functions
  static IDDC T(IDDContent* g, IDDContent* p, IDDContent* d);                // res : address of n for n = g + xp*d (global memo function)
  inline static IDDC T(unsigned int n) {return T((IDDContent*)n, NULL, NULL);}
  static IDDC C(const IDDC &g, const IDDC &p, const IDDC &d);                // res : Impose max(g,d) < xp (Constructor)    
  static IDDC X(const IDDC &n);                                              // res : 2^n-1 (global memo function)
  static IDDC TWICE(const IDDC &n, hash<IDDC1> *H2x = NULL);                 // res : n*2 (multiplication by 2)
  static IDDC AM(const IDDC &n, const IDDC &m, hash<IDDC2> *Ha);             // res : special n+m (local memo function)  
  
  // My internal functions
  static IDDC Y(const IDDC &n, hash<IDDC1> *Hx2);                            // res : Xn/2
  static IDDC HALF(const IDDC &n, IDDC *r = NULL, hash<IDDC1> *Hx2 = NULL);  // res : n/2 (division by 2)   

  // Operators on IDDC defined in the article
  static int  COMP(const IDDC &a, const IDDC &b);                            // res : 1 if a>b, 0 if a==b, -1 if a<b (Comparaison)  
  static IDDC DEC (const IDDC &n);                                           // res : n-1 (Decrement)
  static IDDC INC (const IDDC &n);                                           // res : n+1 (Increment)
  static IDDC RMSB(const IDDC &n, IDDC *i = NULL);                           // res : n-2^i with res<2^i  (WARNING: n != 0) (Remove the MSB)
  static IDDC AMSB(const IDDC &n, const IDDC &i);                            // res : n+2^i (WARNING n<2^i) (Add the MSB)
  static IDDC A   (const IDDC &a, const IDDC &b, hash<IDDC2> *Ha = NULL);    // res : a+b (Addition)
  static IDDC S   (const IDDC &a, const IDDC &b);                            // res : a-b (WARNING: a >= b) (Substraction)
  static IDDC P   (const IDDC &a, const IDDC &b, hash<IDDC2> *Hp = NULL);    // res : a*b (Product)

  // My Operators
  static IDDC NEG (const IDDC &n, const IDDC &q);                            // res : two's complement of n to Xq (WARNING: n<Xq) 
  static IDDC LSH (const IDDC &n, const IDDC &m);                            // res : n << m (Left shift)  
  static IDDC RSH (const IDDC &n, const IDDC &m);                            // res : n >> m (Riht shift)
  static IDDC POW (const IDDC &n, const IDDC &m, hash<IDDC2> *Hp = NULL);    // res : n^m (Power)    
  static inline bool IS_EVEN(IDDContent *n) {if(n->ul_fit()) return n->val()%2==0; else return IS_EVEN(n->g);}
  static int SIZE(IDDContent *n, hash<IDDC1> *visit = NULL);
  
  // IDD conversion
  static inline IDDC OF(unsigned int n) {return T(n);}
  static unsigned int TO(const IDDC &n, int l, bool dec);
  static IDDC OF_MPZ(const mpz_t &n);
  static void TO_MPZ(mpz_t &rop, const IDDC &n);    
 
  public:
  IDDContent(IDDContent* _g, IDDContent* _p, IDDContent* _d);
  ~IDDContent();
  void print(bool rec = false, bool n_only = false);
  bool is_corrupted(mpz_t *n = NULL);
  inline double bitsize() {int s = SIZE(this); return (double)s*(ceil(log2(s+1.))-1.);}
  inline double density() {IDDC l; RMSB(this,&l); return bitsize()/(double)(TO(l,32,1));}  
  
  // Friends
  friend class IDDC;
  
  friend class uIDD;  
  friend bool operator == (const uIDD &n, const unsigned int &m);
  friend bool operator == (const uIDD &n, const uIDD &m);
  friend bool operator != (const uIDD &n, const unsigned int &m);
  friend bool operator != (const uIDD &n, const uIDD &m);
  friend bool operator <  (const uIDD &n, const uIDD &m);
  friend bool operator >  (const uIDD &n, const uIDD &m);
  friend bool operator <= (const uIDD &n, const uIDD &m);
  friend bool operator >= (const uIDD &n, const uIDD &m);
  friend uIDD operator +  (const uIDD &n, const uIDD &m);
  friend uIDD operator -  (const uIDD &n, const uIDD &m);
  friend uIDD operator *  (const uIDD &n, const uIDD &m);
  friend uIDD operator *  (const uIDD &n, const unsigned int &m);
  friend uIDD operator *  (const unsigned int &n, const uIDD &m);   
  friend uIDD operator /  (const uIDD &n, const uIDD &m);
  friend uIDD operator << (const uIDD &n, const uIDD &k);
  friend uIDD operator >> (const uIDD &n, const uIDD &k); 
  
  friend int idd_sign(const uIDD& n);
  friend int idd_compare(const uIDD &n, const uIDD &m); 
  friend uIDD idd_neg_2(const uIDD &n, const uIDD &q);
  friend uIDD idd_pow(const uIDD &n, const uIDD& m);
  friend unsigned int idd_to_int(const uIDD &n);  
  friend void idd_print(const uIDD &n);
    
  friend class sIDD;  
  friend int idd_sign(const sIDD& n);     
  friend int idd_compare(const sIDD& n, const sIDD& m);
  friend sIDD idd_pow(const sIDD &n, const uIDD& m);
  friend int idd_to_int(const sIDD &n);
  friend void idd_print(const sIDD &n);    
  
  // For testing
  friend class IDDTest;  
  friend int main (int argc, char *argv[]);
};



//----------------------------------------------------------------------------//
//                                                                            //
//                               Unsigned IDDC                                //
//                                                                            //
//----------------------------------------------------------------------------//   
class uIDD
{
  private:
  IDDC a;
     
  public:  
  inline uIDD() {a = IDDContent::zero;}
  inline uIDD(unsigned int n) {a = IDDContent::OF(n);}  
  inline uIDD(IDDC _a) {a = _a;}    
  uIDD(const mpz_t &n) {mpz_t np; mpz_init(np); mpz_abs(np,n); a = IDDContent::OF_MPZ(np);}  
  
  // Friend functions & classes
  friend class IDDTest;    
  friend int main (int argc, char *argv[]);
  friend int  idd_sign(const uIDD& n);  
  friend int  idd_compare(const uIDD &n, const uIDD &m);  
  friend uIDD idd_neg_2(const uIDD &n, const uIDD &q);    // negation (two's comp#include<stdlib.h>lement) of n to Xq.
  friend uIDD idd_pow(const uIDD &n, const uIDD& m);  
  friend sIDD idd_pow(const sIDD &n, const uIDD& m);
  friend void idd_print(const uIDD &n);  
   
  // Conversion
  friend unsigned int idd_to_int(const uIDD &n);
  inline void to_mpz(mpz_t &rop) const {IDDContent::TO_MPZ(rop,a);}
  
  // Affectation
  uIDD& operator = (const uIDD &n) {a = n.a; return *this;}  
  uIDD& operator = (const unsigned int &n) {a = IDDContent::OF((unsigned int)labs(n)); return *this;}
  uIDD& operator = (const mpz_t &n) {mpz_t np; mpz_init(np); mpz_abs(np,n); a = IDDContent::OF_MPZ(n); return *this;}

  // Comparison operators     
  friend bool operator == (const uIDD &n, const unsigned int &m);
  friend bool operator == (const uIDD &n, const uIDD &m);
  friend bool operator != (const uIDD &n, const unsigned int &m);  
  friend bool operator != (const uIDD &n, const uIDD &m);  
  friend bool operator <  (const uIDD &n, const uIDD &m);
  friend bool operator >  (const uIDD &n, const uIDD &m);
  friend bool operator <= (const uIDD &n, const uIDD &m);
  friend bool operator >= (const uIDD &n, const uIDD &m);
  
  // Arithmetic operators
  inline void operator ++ () {a = IDDContent::INC(a);}
  inline void operator -- () {a = IDDContent::DEC(a);}
  friend uIDD operator + (const uIDD &n, const uIDD &m);
  friend uIDD operator - (const uIDD &n, const uIDD &m);
  friend uIDD operator * (const uIDD &n, const uIDD &m);
  friend uIDD operator * (const uIDD &n, const unsigned int &m);
  friend uIDD operator * (const unsigned int &n, const uIDD &m);  
  friend uIDD operator / (const uIDD &n, const uIDD &m);
  
  // Logical operators
  friend uIDD operator << (const uIDD &n, const uIDD &k);
  friend uIDD operator >> (const uIDD &n, const uIDD &k);
};
//------------------------------------------------------------------------------
inline bool operator == (const uIDD &n, const unsigned int &m) {return IDDContent::COMP(n.a, IDDContent::OF(m)) == 0;}
inline bool operator == (const uIDD &n, const uIDD &m) {return IDDContent::COMP(n.a, m.a) == 0;}
inline bool operator != (const uIDD &n, const unsigned int &m) {return IDDContent::COMP(n.a, IDDContent::OF(m)) != 0;}
inline bool operator != (const uIDD &n, const uIDD &m) {return IDDContent::COMP(n.a, m.a) != 0;}
inline bool operator <  (const uIDD &n, const uIDD &m) {return IDDContent::COMP(n.a, m.a) < 0;}
inline bool operator >  (const uIDD &n, const uIDD &m) {return IDDContent::COMP(n.a, m.a) > 0;} 
inline bool operator <= (const uIDD &n, const uIDD &m) {return IDDContent::COMP(n.a, m.a) <= 0;}
inline bool operator >= (const uIDD &n, const uIDD &m) {return IDDContent::COMP(n.a, m.a) >= 0;}  
inline uIDD operator +  (const uIDD &n, const uIDD &m) {return uIDD(IDDContent::A(n.a, m.a));}
inline uIDD operator -  (const uIDD &n, const uIDD &m) {return uIDD(IDDContent::S(n.a, m.a));}
inline uIDD operator *  (const uIDD &n, const uIDD &m) {return uIDD(IDDContent::P(n.a, m.a));}
inline uIDD operator *  (const uIDD &n, const unsigned int &m) {return uIDD(IDDContent::P(n.a, IDDContent::OF(m)));}
inline uIDD operator *  (const unsigned int &n, const uIDD &m) {return uIDD(IDDContent::P(IDDContent::OF(n), m.a));}
inline uIDD operator << (const uIDD &n, const uIDD &k) {return uIDD(IDDContent::LSH(n.a, k.a));}
inline uIDD operator >> (const uIDD &n, const uIDD &k) {return uIDD(IDDContent::RSH(n.a, k.a));}


//----------------------------------------------------------------------------//
//                                                                            //
//                                Signed IDDC                                 //
//                                                                            //
//----------------------------------------------------------------------------//
class sIDD   
{
  private:
  bool sign;         // 1 if z<0 else 0  
  IDDC a;
    
  public:
  inline sIDD() {sign = false; a = IDDContent::zero;}
  inline sIDD(const unsigned int &n) {sign = false; a = IDDContent::OF(n);}  
  inline sIDD(const signed int &n) {sign = n<0; a = IDDContent::OF((unsigned int)labs(n));}  
  inline sIDD(unsigned int n) {sign = false; a = IDDContent::OF((unsigned int)n);}  
  inline sIDD(signed int n) {sign = n<0; a = IDDContent::OF((unsigned int)labs(n));}    
  inline sIDD(bool _sign, IDDC _a) {sign = _sign; a = _a;}
  sIDD(const mpz_t &n) {mpz_t np; mpz_init(np); mpz_abs(np,n); sign = mpz_sgn(n)<0; a = IDDContent::OF_MPZ(np);}  

  // Friend functions
  friend class IDDTest;    
  friend int idd_sign(const sIDD &n);
  friend int idd_compare(const sIDD &n, const sIDD& m);
  friend sIDD idd_pow(const sIDD &n, const uIDD& m);
  friend void idd_print(const sIDD &n);
      
  // Conversion
  friend int idd_to_int(const sIDD &n);
  inline void to_mpz(mpz_t &rop) {if(sign) {IDDContent::TO_MPZ(rop,a); mpz_neg(rop, rop);} else IDDContent::TO_MPZ(rop,a);}
  
  // Affectation
  inline void operator = (const sIDD &n) {sign = n.sign; a = n.a;}
  inline void operator = (const signed int &n) {sign = n<0; a = IDDContent::OF((unsigned int)labs(n));}
  void operator = (const mpz_t &n) {mpz_t np; mpz_init(np); mpz_abs(np,n); sign = mpz_sgn(n)<0; a = IDDContent::OF_MPZ(n);}

  // Comparison operators
  inline bool operator == (const int &m)  const {return sign==(m<0) && IDDContent::COMP(a, IDDContent::OF(m)) == 0;}  
  inline bool operator == (const sIDD &m) const {return sign==m.sign && IDDContent::COMP(a, m.a) == 0;}
  inline bool operator != (const int &m)  const {return sign!=(m<0) || IDDContent::COMP(a, IDDContent::OF(m)) != 0;}  
  inline bool operator != (const sIDD &m) const {return sign!=m.sign || IDDContent::COMP(a, m.a) != 0;}
  inline bool operator <  (const sIDD &m) const {int c = IDDContent::COMP(a, m.a); return ((sign&&m.sign)?-c:((sign||m.sign)?(sign?-1:1):c))< 0;}
  inline bool operator >  (const sIDD &m) const {int c = IDDContent::COMP(a, m.a); return ((sign&&m.sign)?-c:((sign||m.sign)?(sign?-1:1):c))> 0;}
  inline bool operator <= (const sIDD &m) const {int c = IDDContent::COMP(a, m.a); return ((sign&&m.sign)?-c:((sign||m.sign)?(sign?-1:1):c))<=0;}
  inline bool operator >= (const sIDD &m) const {int c = IDDContent::COMP(a, m.a); return ((sign&&m.sign)?-c:((sign||m.sign)?(sign?-1:1):c))>=0;}
  
  // Arithmetic operators
  void operator ++ ();
  void operator -- ();
  sIDD operator + (const sIDD &m) const;
  sIDD operator - (const sIDD &m) const;  
  inline sIDD operator * (const sIDD &m) const {return sIDD(sign^m.sign, IDDContent::P(a, m.a));}
};
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
inline int idd_sign(const uIDD &n)
{ 
  return n.a == IDDContent::zero ? 0 : 1;
}
//------------------------------------------------------------------------------
inline int idd_sign(const sIDD &n)
{ 
  return n.a == IDDContent::zero ? 0 : n.sign ? -1 : 1;
}
//------------------------------------------------------------------------------
inline uIDD idd_pow(const uIDD &n, const uIDD& m)
{
  return uIDD(IDDContent::POW(n.a,m.a));
}
//------------------------------------------------------------------------------
inline sIDD idd_pow(const sIDD &n, const uIDD& m)
{
  if(!n.sign)
    return sIDD(false, IDDContent::POW(n.a,m.a));
  else
    return sIDD(!IDDContent::IS_EVEN(*(m.a)), IDDContent::POW(n.a,m.a));
}
//------------------------------------------------------------------------------
inline unsigned int idd_to_int(const uIDD &n)
{
  return IDDContent::TO(n.a,sizeof(int)*8,true);
}
//------------------------------------------------------------------------------
inline int idd_to_int(const sIDD &n)
{
  return IDDContent::TO(n.a,sizeof(int)*8,!n.sign)*(n.sign?-1:1);
}
//------------------------------------------------------------------------------
inline void idd_print(const uIDD &n)
{
  n.a.a->print(false, true);
}
//------------------------------------------------------------------------------
inline void idd_print(const sIDD &n)
{
  if(n.sign) printf("-");
  n.a.a->print(false, true);  
}
//------------------------------------------------------------------------------



//----------------------------------------------------------------------------//
//                                                                            //
//                                 HashTables                                 //
//                                                                            //
//----------------------------------------------------------------------------//
template <class key, class data>
IDDContent::hash<key,data>::hash(unsigned int size, bool (*resizing)(), unsigned int max)
{
  onResizing = resizing;
  if(max == 0) max_size = 0; else max_size = next_prime(max);
  size = next_prime(size);
  if(max_size!=0 && size>max_size) size = max_size;
  init_size = size;
  T.resize(size);
  n_items = 0;
}
//------------------------------------------------------------------------------
template <class key, class data>
IDDContent::hash<key, data>::~hash()
{
  for(unsigned int i=0; i<T.size(); i++)
    if(T[i]) {delete T[i];}
}
//------------------------------------------------------------------------------
template <class key, class data>
int IDDContent::hash<key, data>::next_prime(int n)
{
  if(n%2 == 0) return next_prime(n+1);
  int sq = sqrt(n);
  int i = 3;
  for(;i<=sq; i+=2)
    if(n%i == 0)
      break;
  if(i>sq)
    return n;
  else
    return next_prime(n+2);
}
//------------------------------------------------------------------------------
template <class key, class data>
bool IDDContent::hash<key, data>::exists(const key &k)
{
  int i = k.hash() % T.size(); 
  if(!T[i])
    return false;
  else
    return T[i]->find(k) != T[i]->end();
}
//------------------------------------------------------------------------------
template <class key, class data>
data IDDContent::hash<key,data>::find(const key &k)
{
  int i = k.hash() % T.size(); 
  return (*T[i])[k];
}
//------------------------------------------------------------------------------
template <class key, class data>
void IDDContent::hash<key,data>::insert(const key &k, const data &d)
{ 
  if(T.size()!=max_size && n_items >= T.size())
    resize(T.size()*2);      
  if(n_items>=T.size())   //full
    return;  
  
  n_items++;
  int i = k.hash() % T.size();  
  
  if(!T[i])
  {
    T[i] = new map<key, data>;
  }
  (*T[i])[k] = d;  
}
//------------------------------------------------------------------------------
template <class key, class data>
void IDDContent::hash<key,data>::resize(unsigned int size)
{
  if(onResizing) 
  {
      if(onResizing()) return;
    else
      if(clean()) return;
  }
  if(max_size!=0 && size>max_size)
  {
    if(T.size()==max_size)
      return;
    else
      size = max_size;
  }
  else
    size = next_prime(size);  
  vector<map<key, data>* > S;
  S.resize(size);
  
  for(unsigned int i=0; i<T.size(); i++)
    if(T[i])
    {
      typename map<key, data>::iterator it = T[i]->begin();
      while(it != T[i]->end())
      {
        int j = it->first.hash() % size;
        if(!S[j])
        {
          S[j] = new map<key, data>;
        }
        (*S[j])[it->first] = it->second;
        it++;
      }
      delete T[i];
    }
  T = S;
}
//------------------------------------------------------------------------------
template <class key, class data>
void IDDContent::hash<key,data>::remove(const key &k)
{
  int i = k.hash() % T.size();
  if(T[i])
  {
    typename map<key, data>::iterator it = T[i]->find(k);
    
    if(it != T[i]->end())
    {
      n_items--;  
      T[i]->erase(it);
    }
  }
}
//------------------------------------------------------------------------------
template <class key, class data>
bool IDDContent::hash<key,data>::clean()
{   
  for(unsigned int i=0; i<T.size();i++)
    if(T[i])
    {
      typename map<key,data>::iterator it;
      list<key> rm;
      for(it = T[i]->begin(); it!=T[i]->end(); ++it)
        if(it->second->n_shared==0) 
          rm.push_back(it->first);
      typename list<key>::iterator iter;
      for(iter=rm.begin(); iter!=rm.end(); ++iter)
        if(T[i]->find(*iter) != T[i]->end())
          T[i]->erase(*iter);
    }  

  return items_count() < size()/2;
}
//------------------------------------------------------------------------------
template <class key, class data>
void IDDContent::hash<key,data>::clear()
{
  for(unsigned int i=0; i<T.size(); i++)
  {
    if(T[i]) 
    {
      delete T[i]; 
    }
    T[i] = NULL;   
  }
  T.resize(init_size);
}
//------------------------------------------------------------------------------
template <class key, class data>
void IDDContent::hash<key,data>::show_stats()
{
  double var=0, ect;
  double mean = (double)n_items/(double)T.size();
  double max_diff=-1;
  int max_diff_index=-1;
  unsigned int mem = sizeof(T) + T.size()*sizeof(map<key,data>*);
  for(unsigned int i=0;i<T.size();i++)
  {
    double diff;
    if(T[i])
    {
      diff = fabs((double)T[i]->size()-mean);
      mem += sizeof(map<key,data>)+T[i]->size()*(sizeof(key)+sizeof(data));
    }
    else
      diff = mean;
    if(diff>max_diff)
    {
      max_diff=diff;
      max_diff_index = i;
    }  
    var += pow(diff,2);
  }
  var /= T.size();
  ect = sqrt(var);
  
  printf("Table items: %ld/%d\n", n_items, T.size());
  printf("Mean bucket size: %f\n", mean);
  printf("Standard deviation: %f\n", ect); 
  printf("Max diff on index %d: %f\n",max_diff_index, max_diff);
  printf("Memory estimate: %ld\n",mem);  
}
//------------------------------------------------------------------------------
template <class key, class data>
unsigned int IDDContent::hash<key,data>::memory()
{
  unsigned int mem = sizeof(T) + T.size()*sizeof(map<key,data>*);
  for(unsigned int i=0;i<T.size();i++)
    if(T[i])
      mem += sizeof(map<key,data>)+T[i]->size()*(sizeof(key)+sizeof(data));
  return mem;  
}
//------------------------------------------------------------------------------
#endif
