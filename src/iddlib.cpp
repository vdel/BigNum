#include<assert.h>
#include <limits.h>

#include "iddlib.h"

#include "monitor.h"
extern int call_constructor, call_clean_Ht;
extern int call_T, call_int_T;
extern int call_COMP, call_C;
extern int call_OF, saved_call_OF;
extern int call_DEC, call_INC, call_int_INC;
extern int call_X, saved_call_X;
extern int call_RMSB, call_AMSB;
extern int call_TWICE, saved_call_TWICE;
extern int call_A, saved_call_A, call_int_A;
extern int call_P, saved_call_P, call_int_P;
extern int call_NEG;
extern int call_LSH, call_RSH;
extern int call_POW, call_Y, call_HALF;

IDDContent::hash<IDDContent::IDDC3>  IDDContent::Ht = IDDContent::hash<IDDContent::IDDC3>(66161, clean_Ht);
IDDContent::hash<IDDContent::IDDC1>  IDDContent::Hx = IDDContent::hash<IDDContent::IDDC1>(1024);
IDDC IDDContent::zero = IDDC(new IDDContent((IDDContent*)0, NULL, NULL));
IDDC IDDContent::one  = IDDC(new IDDContent((IDDContent*)1, NULL, NULL));
IDDC IDDContent::Xws  = IDDC(new IDDContent((IDDContent*)IDD_WORD_SIZE, NULL, NULL));


//----------------------------------------------------------------------------//
//                                                                            //
//                         Protected Container                                //
//                                                                            //
//----------------------------------------------------------------------------//
IDDC::IDDC() {a = IDDContent::zero.a; a->n_shared++;}
IDDC::IDDC(const IDDC &s) {a = s.a; a->n_shared++;}
IDDC::IDDC(IDDContent* _a) {a=_a; a->n_shared++;}
IDDC::~IDDC() {a->n_shared--;}
IDDC& IDDC::operator = (const IDDC &m)
{
  m.a->n_shared++;
  a->n_shared--;  
  a = m.a;
  return *this;
}
IDDC& IDDC::operator = (IDDContent *n)
{
  n->n_shared++;
  a->n_shared--;  
  a = n;
  return *this;
}



//----------------------------------------------------------------------------//
//                                                                            //
//                               IDD Container                                //
//                                                                            //
//----------------------------------------------------------------------------//
IDDContent::IDDContent(IDDContent* _g, IDDContent* _p, IDDContent* _d):g(_g),p(_p),d(_d),n_shared(0)
{
  #ifdef MONITORING
  call_constructor++;
  #endif
  if(!ul_fit())
  {
    g->n_shared++;
    p->n_shared++;
    d->n_shared++;
  }
}
//------------------------------------------------------------------------------
IDDContent::~IDDContent()
{
  Ht.remove(IDDC3(g,p,d));
   
  if(!ul_fit())
  {
    g->n_shared--; if(!g->n_shared) delete g;
    p->n_shared--; if(!p->n_shared) delete p;
    d->n_shared--; if(!d->n_shared) delete d;
  }
}
//------------------------------------------------------------------------------
bool IDDContent::clean_Ht()
{   
  #ifdef MONITORING
  call_clean_Ht++;
  #endif
  
  for(unsigned int i=0; i<Ht.T.size();i++)
    if(Ht.T[i])
    {
      map<IDDC3,IDDContent*>::iterator it;
      list<pair<IDDC3,IDDContent*> > rm;
      for(it = Ht.T[i]->begin(); it!=Ht.T[i]->end(); ++it)
      {
        if(it->second->n_shared==0) 
          rm.push_back(pair<IDDC3,IDDContent*>(it->first, it->second));
      }
      list<pair<IDDC3,IDDContent*> >::iterator iter;
      for(iter=rm.begin(); iter!=rm.end(); ++iter)
        if(Ht.T[i]->find(iter->first) != Ht.T[i]->end())
        {
          delete iter->second;
        }
    }  

  return Ht.items_count() < Ht.size()/2;
}
//------------------------------------------------------------------------------
unsigned int IDDContent::TO(const IDDC &n, int l, bool dec)
{
       if(n == zero) return 0;
  else if(n == one) return 1;
  else if(n->ul_fit()) return n->val();
  else
  {
    if(COMP(AMSB(zero,OF((unsigned int)l)), n) > (dec?0:-1))
    {
      unsigned int g = TO(n->g, l, dec);
      unsigned int p = TO(n->p, l, dec);
      unsigned int d = TO(n->d, l, dec);
      p = 1 << (1 << p);
      return g + p*d;
    }  
    else
      return TO(DEC(AMSB(zero,OF((unsigned int)l))), l, dec);
  }
}
//------------------------------------------------------------------------------
IDDC IDDContent::OF_MPZ(const mpz_t &n)
{
  IDDC res = zero;
  unsigned int i = -1;
  while((i=mpz_scan1(n, i+1)) != ULONG_MAX)
    res = AMSB(res, OF(i));
  return res;
}
//------------------------------------------------------------------------------
void IDDContent::TO_MPZ(mpz_t &rop, const IDDC &n)
{
       if(n == zero) mpz_set_ui(rop, 0);
  else if(n == one) mpz_set_ui(rop, 1);
  else
  {
    if(n->ul_fit())
      mpz_set_ui(rop,n->val());
    else    
    {
      mpz_t _1;
      mpz_t g,p,d;
      mpz_init_set_ui(_1, 1);
      mpz_init(g);
      mpz_init(p);
      mpz_init(d);
      TO_MPZ(g,n->g);
      unsigned int pui = TO(n->p, 32, 0);
      TO_MPZ(d,n->d);
      pui = 1<<pui;
      mpz_mul_2exp(p, _1, pui);
      mpz_mul(rop,p,d);
      mpz_add(rop,g,rop);
    }
  }
}
//------------------------------------------------------------------------------
IDDC IDDContent::T(IDDContent* g, IDDContent* p, IDDContent* d)
{ 
  #ifdef MONITORING
  if(p) call_T++; else call_int_T++;
  #endif 

  if(p && p->ul_fit() && p->val()<IDD_WORD_SIZE)
  {
    g = (IDDContent*)(g->val()+(1<<(1<<p->val()))*d->val());
    p = NULL;
    d = NULL;
  }
  if(!p)
  {
    if(g==(IDDContent*)0) return zero;
    if(g==(IDDContent*)1) return one;
    if(g==(IDDContent*)IDD_WORD_SIZE) return Xws;
  }
  
  IDDC3 key(g, p, d);
  if(Ht.exists(key))
    return Ht.find(key);
  else
  {
    IDDC n(new IDDContent(g, p, d));
    Ht.insert(key, *n);
    return n;
  }
}
//------------------------------------------------------------------------------
IDDC IDDContent::C(const IDDC &g, const IDDC &p, const IDDC &d)
{
  #ifdef MONITORING
  call_C++;
  #endif 
  if(d == zero) return g;   // d == 0 ?

  // p == d->p ??
  if(!d->ul_fit())   
  {
    if(p == d->p) return C(C(g, p, d->g), INC(p), d->d); // p == d.p
  }
  else if(d != one && p->ul_fit() && p->val()<IDD_WORD_SIZE)
  {
    unsigned int dp = d->get_p();
    if(p->val() == dp)
      return C(C(g, p, T(d->get_g(dp))), INC(p), T(d->get_d(dp))); // p == d.p
  }
  // p != d->p
  
  // p == d.g ??
  if(!g->ul_fit())
  {
    if(p == g->p) return C(g->g, p, A(d, g->d)); // p == d.g
  }
  else if(g != zero && g != one && p->ul_fit() && p->val()<IDD_WORD_SIZE)
  {
    unsigned int gp = g->get_p();
    if(p->val() == gp)   
      return C(T(g->get_g(gp)), p, A(d, T(g->get_d(gp)))); // p == d.g
  }
  // p != g->p
  
  return T(*g, *p, *d);
}
//------------------------------------------------------------------------------
int IDDContent::COMP(const IDDC &a, const IDDC &b)
{
  #ifdef MONITORING
  call_COMP++;
  #endif
  if(a->ul_fit())
  {
    if(!b->ul_fit())
      return -1;
    else
      return a->val()==b->val()?0:(a->val()>b->val()?1:-1);
  }
  else if(b->ul_fit())
    return 1;
  else
    return a == b       ? 0 :             
           a->p != b->p ? COMP(a->p, b->p) :
           a->d != b->d ? COMP(a->d, b->d) :
                          COMP(a->g, b->g);
}
//------------------------------------------------------------------------------
IDDC IDDContent::DEC(const IDDC &n)
{
  #ifdef MONITORING
  call_DEC++;
  #endif
  if(n->ul_fit())
  {
    if(n->val())
      return T(n->val()-1);
    else 
      return n;
  }
  else
    return n->g != zero ? T(*DEC(n->g), n->p, n->d) :
           n->d == one  ? X(n->p) :
                          T(*X(n->p), n->p, *DEC(n->d));
}
//------------------------------------------------------------------------------
IDDC IDDContent::INC(const IDDC &n)
{
  if(n->ul_fit())
  {
    #ifdef MONITORING
    call_int_INC++;
    #endif  
    if(n->val()!=ULONG_MAX)
      return T(n->val()+1);
    else 
      return T(*zero, *Xws, *one);
  }
  else  
  {
    #ifdef MONITORING
    call_INC++;
    #endif
    return n->g != X(n->p) ? T(*INC(n->g), n->p, n->d) :
           n->d != X(n->p) ? T(*zero, n->p, *INC(n->d)) :
                             T(*zero, *INC(n->p), *one);
  }
}
//------------------------------------------------------------------------------
IDDC IDDContent::X(const IDDC &n)
{
  IDDC1 key(*n);
  IDDC res;
  if(Hx.exists(key))
  {
    res = Hx.find(key);
    #ifdef MONITORING
    saved_call_X++;
    #endif    
  }
  else
  {
    #ifdef MONITORING
    call_X++;
    #endif  
    if(n->ul_fit() && n->val()<=IDD_WORD_SIZE)
    {
      if(n->val() == IDD_WORD_SIZE)
        res = T(-1);
      else
        res = T((1<<(1<<n->val()))-1);
    }
    else
    {
      IDDC q = DEC(n);
      IDDC rec = X(q);
      res = T(*rec, *q, *rec);
    }
    res->n_shared++;
    Hx.insert(key, *res);
  }
  return res;  
}
//------------------------------------------------------------------------------
IDDC IDDContent::RMSB(const IDDC &n, IDDC *i)
{
  if(n->ul_fit())
  {
    int k = -1;
    unsigned int m = n->val();
    assert(m!=0);
    while(m)
    {
      k++;
      m=m>>1;
    }
    if(i) *i = T(k);
    return T(n->val()-(1<<k));
  }
  else
  {
    #ifdef MONITORING  
    call_RMSB++;
    #endif
    IDDC l;
    IDDC e = RMSB(n->d, &l);
    if(i) *i = AMSB(l, n->p);
    return C(n->g, n->p, e);
  }
}
//------------------------------------------------------------------------------
IDDC IDDContent::AMSB(const IDDC &n, const IDDC &i)
{ 
  if(i->ul_fit() && i->val()<(1<<IDD_WORD_SIZE))
  {
    unsigned int k = (unsigned int)1<<i->val();
    assert(n->ul_fit() && n->val()<k);
    return T(n->val()+k);   
  }
  else
  {
    #ifdef MONITORING
    call_AMSB++;
    #endif
    IDDC l;
    IDDC e = RMSB(i, &l);
    unsigned int np_ul = n->get_p();
    IDDC np = n->ul_fit() ? T(n->get_p()) : n->p;
    if(n == zero || n == one || COMP(l, np) > 0)
      return T(*n, *l, *AMSB(zero, e));
    else
    {
      IDDC ng = n->ul_fit() ? T(n->get_g(np_ul)) : n->g;
      IDDC nd = n->ul_fit() ? T(n->get_d(np_ul)) : n->d;      
      return T(*ng, *l, *AMSB(nd, e));
    }
  }
}
//------------------------------------------------------------------------------
IDDC IDDContent::TWICE(const IDDC &n, hash<IDDC1> *H2x)
{
  if(n->ul_fit())
  {
    if(n->val()&0x80000000)
    {
      IDDC tmp = T(n->val()<<1);
      return T(*tmp,*Xws,*one);
    }
    else
      return T(n->val()<<1);
  }
  else  
  {
    bool alloc = false;
    if(H2x == NULL)
    {
      H2x = new hash<IDDC1>(1031);
      alloc = true;
    }
    
    IDDC1 key(*n); 
    IDDC res;
    
    if(H2x->exists(key))
    {
      res = H2x->find(key);
      #ifdef MONITORING
      saved_call_TWICE++;
      #endif
    }
    else
    {
      #ifdef MONITORING
      call_TWICE++;
      #endif
      res = C(TWICE(n->g, H2x), n->p, TWICE(n->d, H2x));
      H2x->insert(key, *res);
    }

    if(alloc)
    {
      delete H2x;
    }
    return res;  
  }
}
//------------------------------------------------------------------------------
IDDC IDDContent::A(const IDDC &a, const IDDC &b, hash<IDDC2> *Ha)
{
  if(a->ul_fit() && b->ul_fit())
  {
    #ifdef MONITORING    
    call_int_A++;     
    #endif  
    if(a->val()<=ULONG_MAX-b->val()) 
      return T(a->val()+b->val());
    else  // a+b = xM + (a+b-xM)
    {
      IDDC tmp = T(a->val()-(ULONG_MAX-b->val())-1);
      return T(*tmp,*Xws,*one);
    }
  }
  else if(COMP(a, b) > 0)
    return A(b,a,Ha);
  else
  {
    #ifdef MONITORING    
    call_A++;     
    #endif

    bool alloc = false;
    if(Ha == NULL)
    {
    
      Ha = new hash<IDDC2>(1031);
      alloc = true;
    }
    
    IDDC res;
    
    if(a == b)
      res = TWICE(a);
    else     
    // if b fits, then a fits and we cannot be here, so b does not fit.
    if(a->ul_fit())  // then a->p<b->p
      res = C(A(a, b->g, Ha), b->p, b->d);
    else
      res = COMP(a->p, b->p) < 0 ? C(A(a, b->g, Ha), b->p, b->d) :
                                   AM(a, b, Ha);    

    if(alloc)
    {
      delete Ha;
    }
      
    return res;
  }
}
//------------------------------------------------------------------------------
IDDC IDDContent::AM(const IDDC &a, const IDDC &b, hash<IDDC2> *Ha)
{
  // By construction, neither a nor b fits.
  IDDC2 key(*a, *b);
  IDDC res;
  
  if(Ha->exists(key))
  {
    res = Ha->find(key);
    #ifdef MONITORING    
    saved_call_A++;
    #endif
  }
  else
  {
    res = C(A(a->g, b->g, Ha), a->p, A(a->d, b->d, Ha));
    Ha->insert(key, *res);     
  }

  return res;  
}
//------------------------------------------------------------------------------
IDDC IDDContent::S(const IDDC &a, const IDDC &b)
{
  assert(COMP(a, b) >= 0);
  if(a->ul_fit())
    return T(a->val()-b->val());
  else
    return b == zero ? a :
           b == one  ? DEC(a) :
           b == a    ? zero :
                       A(INC(a), NEG(b, INC(a->p)) )->g;
}
//------------------------------------------------------------------------------
IDDC IDDContent::P(const IDDC &a, const IDDC &b, hash<IDDC2> *Hp)
{    
  if(COMP(a,b) > 0)
    return P(b, a, Hp);
  else 
  { 
    bool alloc = false;
    if(Hp == NULL)
    {    
      Hp = new hash<IDDC2>(1031, NULL, 66161);
      alloc = true;
    } 
  
    IDDC res;
    IDDC2 key(*a, *b);

    if(Hp->exists(key))
    {
      res = Hp->find(key);      
      #ifdef MONITORING      
      saved_call_P++;
      #endif
    }
    else
    {   
      if(a->ul_fit() && b->ul_fit())
      {
        #ifdef MONITORING    
        call_int_P++;
        #endif
        uint64_t mul = uint64_t(a->val())*uint64_t(b->val());
        unsigned int d = mul>>32;
        res = T(mul&0xffffffff);
        if(d)
        {
          IDDC tmp = T(d);
          res = T(*res,*Xws,*tmp); 
        }              
      }
      else
      {
        #ifdef MONITORING    
        call_P++;
        #endif        
        // b does not fit because otherwise a fits too and we would not be here.
        IDDC g = P(a, b->g, Hp);
        IDDC d = P(a, b->d, Hp);
              
        bool overflow_g = !g->ul_fit() && COMP(g->p,b->p) > 0; 
        bool overflow_d = !d->ul_fit() && COMP(d->p,b->p) > 0;

        if(overflow_g && overflow_d)
          res = C(C(g->g,b->p,d->g),d->p,C(g->d,b->p,d->d));
        else if(overflow_g) 
          res = C(C(g->g,b->p,d),g->p,g->d);
        else if(overflow_d) 
          res = C(C(g,b->p,d->g),d->p,C(zero,b->p,d->d));
        else
          res = C(g, b->p, d);   
      }
     
      Hp->insert(key, *res);
    }     

    if(alloc)
    {
      delete Hp;  
    }
    return res;            
  }   
}
//------------------------------------------------------------------------------
IDDC IDDContent::NEG(const IDDC &n, const IDDC &q)
{
  #ifdef MONITORING    
  call_NEG++;
  #endif
  if(n->ul_fit())
  {
    if(q->ul_fit() && q->val()<=IDD_WORD_SIZE)
    {
      unsigned int dec = (1<<IDD_WORD_SIZE) - (1<<q->val());
      unsigned int not_n = ~n->val();
      not_n = (not_n << dec) >> dec;
      return T(not_n);
    }
    else
    {
      IDDC q2 = DEC(q);
      return T(*NEG(n,q2), *q2, *X(q2));
    }
  }
  else
  {
    IDDC q2 = DEC(q);
    if(n->p == q2)
      return C(NEG(n->g, n->p), n->p, NEG(n->d, n->p));
    else
      return T(*NEG(n,q2), *q2, *X(q2));
  }
}
//------------------------------------------------------------------------------
IDDC IDDContent::LSH(const IDDC &n, const IDDC &m)
{
  #ifdef MONITORING    
  call_LSH++;
  #endif
  if(m == zero)
    return n;
  else if(n == zero)
    return zero;
  else if(n == one)
    return AMSB(zero, m);
  else
  {
    IDDC l;  
    IDDC e = RMSB(m, &l);
    unsigned int np_ul = n->get_p();
    IDDC np = n->ul_fit() ? T(np_ul) : n->p;   
    if(COMP(l, np) > 0)
      return C(zero, l, LSH(n, e));
    else
    {
      IDDC ng = n->ul_fit() ? T(n->get_g(np_ul)) : n->g;   
      IDDC nd = n->ul_fit() ? T(n->get_d(np_ul)) : n->d;         
      return C(LSH(ng, m), np, LSH(nd, m));
    }
  }  
}
//------------------------------------------------------------------------------
IDDC IDDContent::RSH(const IDDC &n, const IDDC &m)
{
  #ifdef MONITORING    
  call_RSH++;
  #endif
  if(m == zero)
    return n;
  else if(n->ul_fit())
  {
    if(m->ul_fit() && m->val()<(1<<IDD_WORD_SIZE))
      return T(n->val()>>m->val());
    else
      return zero;    
  }
  else
  {
    IDDC k = AMSB(zero, n->p);
    if(COMP(k, m) >= 0)
    {
      k = S(k, m);
      return A(RSH(n->g, m), LSH(n->d, k));
    }
    else
    {
      k = S(m, k);
      return RSH(n->d, k);
    }
  }
}
//------------------------------------------------------------------------------
/*IDDC IDDContent::D(IDDC a, IDDC b, IDDC *r)
{
  assert(b != zero);
  if(b == one)
  {
    if(r) *r = zero;
    return a;
  }
  else if(COMP(b,a) > 0)
  {
    if(r) *r = a;
    return zero;
  }
  else if(a == b)
  {
    if(r) *r = zero;
    return one;    
  }
  else
  {
    if(COMP(a->d,b) < 0)  // b >= 2
    {
      if(a->p == zero) // a == 2 || a == 3
      {
        // if a==2: b>=2, b!=a -> b >= 3 -> q=0, r=a
        // if a==3: if b==2 then q=1,r=1 else q=0, r=a
        if(COMP(b,a) > 0)
        {
          if(r) *r = a;
          return zero;  
        }
        else
        {
          if(r) *r = one;
          return one;  
        }
      }
      else
      {
        IDDC pm1 = DEC(a->p);
        IDDC rg;
        IDDC qg = D(a->g, b, &rg);
        IDDC rd;
        IDDC qd = D(LSH(a->d, AMSB(zero, pm1)), b, &rd);
        IDDC q  = D(C(rg, pm1, rd), b, r);
        return A(q, C(qg, pm1, qd));
      }
    }
    else
    {
      IDDC rg;
      IDDC qg = D(a->g, b, &rg);
      IDDC rd;
      IDDC qd = D(a->d, b, &rd);
      IDDC q  = D(C(rg, a->p, rd), b, r);
      return A(q, C(qg, a->p, qd));
    }
  }
  return NULL;
}*/
//------------------------------------------------------------------------------
IDDC IDDContent::HALF(const IDDC &n, IDDC *r, hash<IDDC1> *Hx2)
{
  #ifdef MONITORING    
  call_HALF++;
  #endif
  if(n->ul_fit())
  {
    if(r) *r = (n->val()%2) ? one : zero;
    return T(n->val()/2);
  } 
  else
  {
    bool alloc = false;
    if(Hx2 == NULL)
    {
      Hx2 = new hash<IDDC1>(1031);
      alloc = true;
    }

    IDDC1 key(*n);
    IDDC res;
    
    if(Hx2->exists(key))
      res = Hx2->find(key);
    else
    {
      IDDC rd;
      IDDC qd = HALF(IDDC(n->d), &rd, Hx2);    
      IDDC qg = HALF(IDDC(n->g), r, Hx2);
      
      if(rd == zero)
        res = T(*qg, n->p, *qd);
      else
        res = C(A(qg, Y(n->p, Hx2)), n->p, qd);
      Hx2->insert(key, *res);        
    }
          
    if(alloc)
    {
      delete Hx2;
    }
    return res;
  }
}
//------------------------------------------------------------------------------
IDDC IDDContent::Y(const IDDC &n, hash<IDDC1> *Hx2)
{
  #ifdef MONITORING    
  call_Y++;
  #endif
  if(n->ul_fit() && n->val()<IDD_WORD_SIZE)
    return T((1<<(1<<n->val()))/2);
  else
  {
    IDDC1 key(*n);
    IDDC res;
    if(Hx2->exists(key))
      res = Hx2->find(key);
    else
    {
      IDDC prev = DEC(n);
      IDDC tmp = Y(prev, Hx2);
      res = T(*zero, *prev, *tmp);
      Hx2->insert(key, *res);
    }
    return res;  
  }
}
//------------------------------------------------------------------------------
IDDC IDDContent::POW(const IDDC &n, const IDDC &m, hash<IDDC2> *Hp)
{
  #ifdef MONITORING    
  call_POW++;
  #endif
  
  if(m == zero)
    return one;
  else if(m == one)
    return n;
  else
  {
    bool alloc = false;
    if(Hp == NULL)
    {
      Hp = new hash<IDDC2>(1031, NULL, 66161);
      alloc = true;
    }  
  
    IDDC res,r;
    IDDC m2 = HALF(m, &r);
    IDDC n2 = POW(n, m2, Hp);
        
    if(r == zero)
      res = P(n2,n2, Hp);
    else
      res = P(n, P(n2,n2,Hp), Hp);
    
    if(alloc)
    {
      delete Hp;
    }
    return res;
  }
}
//------------------------------------------------------------------------------
void IDDContent::print(bool rec, bool n_only)
{
  mpz_t n; 
  mpz_init(n);
  TO_MPZ(n,this);
  mpz_out_str(stdout, 10, n);
  if(n_only) return;
  
  if(ul_fit())
  {
    cout << " = (fit)";
    if(rec) cout << endl;
  }
  else
  { 
    cout << " = ("; g->print(false,true); cout << ", "; p->print(false,true); cout << ", "; d->print(false,true); cout << ")";
    if(rec)
    {
      cout << endl;
      if(g) g->print(true);
      if(p) p->print(true);  
      if(d) d->print(true);  
    }
  }
}
//------------------------------------------------------------------------------
bool IDDContent::is_corrupted(mpz_t *n)
{
  if(p && ((unsigned int)g<1000 || (unsigned int)p<1000 || (unsigned int)d<1000))
    return true;
  if(ul_fit())
  {
    if(n) mpz_set_ui(*n, val());
    return false;
  }
  else
  {
    mpz_t gval,pval,dval;
    mpz_init(gval);
    mpz_init(pval);
    mpz_init(dval);
    if(g->is_corrupted(&gval) || p->is_corrupted(&pval) || d->is_corrupted(&dval)) 
    {
      mpz_clear(gval);
      mpz_clear(pval);
      mpz_clear(dval);
      return true;
    }
    else
    {
      unsigned int maxui = mpz_get_ui(pval);
      maxui = 1<<maxui;
      mpz_t _1, max;
      mpz_init_set_ui(_1, 1);
      mpz_init(max);
      mpz_mul_2exp(max, _1, maxui);
      if(mpz_cmp(max, gval) <= 0 || mpz_cmp(max, dval) <= 0)
      {
        mpz_clear(gval);
        mpz_clear(pval);
        mpz_clear(dval);
        return true;
      }
      else
      {  
        if(n) 
        {
          mpz_mul(*n,max,dval);
          mpz_add(*n,gval,*n);
        }
        mpz_clear(gval);
        mpz_clear(pval);
        mpz_clear(dval);
        return false;      
      }    
    }
  }
}
//------------------------------------------------------------------------------
int IDDContent::SIZE(IDDContent *n, hash<IDDC1> *visit)
{
  if(n == zero || n == one)
    return 0;
    
  bool alloc = false;
  if(!visit)
  {
    alloc = true;
    visit = new hash<IDDC1>(1024);
  }
 
  int size = 0;
  
  if(!visit->exists(IDDC1(n)))
  {
    visit->insert(IDDC1(n), n);  
    if(n->ul_fit())
    {
      int mpi = n->get_p();
      IDDC mp = T(mpi);
      IDDC mg = T(n->get_g(mp));
      IDDC md = T(n->get_d(mp));
      size = 1 + SIZE(*mg, visit) + SIZE(*mp, visit) + SIZE(*md,visit);
    }  
    else
      size = 1 + SIZE(n->g, visit) + SIZE(n->p, visit) + SIZE(n->d,visit);
  }
    
  if(alloc)
  {
    delete visit;
  }
  return size;
}
//------------------------------------------------------------------------------





//----------------------------------------------------------------------------//
//                                                                            //
//                               Unsigned IDD                                 //
//                                                                            //
//----------------------------------------------------------------------------//   
uIDD idd_neg_2(const uIDD &n, const uIDD &q)
{
  if     (n.a == IDDContent::zero) return uIDD(IDDContent::X(q.a));
  else if(n.a == IDDContent::one)  return uIDD(IDDContent::DEC(IDDContent::X(q.a)));
  else
  {
    assert(IDDContent::COMP(q.a, n.a->p) > 0);
    return uIDD(IDDContent::NEG(n.a, q.a));
  }
}
//------------------------------------------------------------------------------







//----------------------------------------------------------------------------//
//                                                                            //
//                                Signed IDD                                  //
//                                                                            //
//----------------------------------------------------------------------------//-
void sIDD::operator ++ () 
{
  if(sign)
  {
    if(a == IDDContent::one)
    {
      sign = false;
      a = IDDContent::zero;
    }
    else
      a = IDDContent::DEC(a);
  }
  else
    a = IDDContent::INC(a);  
}
//------------------------------------------------------------------------------
void sIDD::operator -- () 
{
  if(!sign)
  {
    if(a == IDDContent::zero)
    {
      sign = true;
      a = IDDContent::one;
    }
    else
      a = IDDContent::DEC(a);
  }
  else
    a = IDDContent::INC(a);  
}
//------------------------------------------------------------------------------
sIDD sIDD::operator + (const sIDD &m) const 
{
  if(sign == m.sign)
    return sIDD(sign, IDDContent::A(a, m.a));
  else if(sign)  // this neg, m pos -> (*m.a) - (*a)
  {
    if(IDDContent::COMP(a, m.a) > 0)
      return sIDD(true, IDDContent::S(a, m.a));
    else
      return sIDD(false, IDDContent::S(m.a, a));
  }
  else  // this pos, m neg -> (*a) - (*m.a)
  {
    if(IDDContent::COMP(a, m.a) >= 0) 
      return sIDD(false, IDDContent::S(a, m.a));
    else
      return sIDD(true, IDDContent::S(m.a, a));    
  }   
}
//------------------------------------------------------------------------------
sIDD sIDD::operator - (const sIDD &m) const 
{
  if(sign != m.sign)
    return sIDD(sign, IDDContent::A(a, m.a));
  else if(sign)  // both neg -> (*m.a) - (*a)
  {
    if(IDDContent::COMP(a, m.a) > 0)
      return sIDD(true, IDDContent::S(a, m.a));
    else
      return sIDD(false, IDDContent::S(m.a, a));
  }
  else  // both pos -> (*a) - (*m.a)
  {
    if(IDDContent::COMP(a, m.a) >= 0)
      return sIDD(false, IDDContent::S(a, m.a));
    else
      return sIDD(true, IDDContent::S(m.a, a));    
  }   
}
//------------------------------------------------------------------------------
int idd_compare(const sIDD& n, const sIDD& m)
{
  if(n.sign != m.sign)
    return n.sign ? -1 : 1;
  else
  {
    int Ncmp = IDDContent::COMP(n.a, m.a);
    return n.sign ? -Ncmp : Ncmp;
  }
}
//------------------------------------------------------------------------------
