!=======================================================================
!  NBODY - integrador de N corpos para satelites em torno de um planeta
!          achatado (J2 e J4), com saida de elementos osculadores e
!          geometricos (Renner & Sicardy 2006).
!
!  Uso:   ./nbody                 (le "simulacao.in")
!         ./nbody outro_arquivo.in
!
!  Os dados dos planetas ficam em "planetas.dat" (ou no arquivo indicado
!  por "arquivo_planetas" no .in).
!
!  Unidades internas (as mesmas do codigo original ips_mas_pra_aei):
!     comprimento = raio equatorial do planeta
!     tempo       = dia
!     "massas"    = G*m  (em R_p^3/dia^2)
!  Referencial: planetocentrico, plano xy = equador do planeta.
!  Integrador: RA15 (Everhart 1985).
!=======================================================================

!-----------------------------------------------------------------------
module constantes
  implicit none
  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi     = 3.141592653589793238_dp
  real(dp), parameter :: pi2    = 2.0_dp*pi
  real(dp), parameter :: deg2rad = pi/180.0_dp
  real(dp), parameter :: rad2deg = 180.0_dp/pi
  real(dp), parameter :: dia    = 86400.0_dp            ! s
  real(dp), parameter :: G_SI_km = 6.67430e-20_dp       ! km^3 kg^-1 s^-2 (CODATA 2018)
end module constantes

!-----------------------------------------------------------------------
module sistema
  ! Estado global do sistema (planeta + satelites)
  use constantes
  implicit none
  ! planeta
  character(len=64) :: nome_planeta = ''
  real(dp) :: gm_planeta_km = 0.0_dp   ! km^3/s^2
  real(dp) :: raio_km       = 0.0_dp   ! km
  real(dp) :: j2 = 0.0_dp, j4 = 0.0_dp
  real(dp) :: gmp                      ! G*M_p em R_p^3/dia^2
  ! satelites
  integer :: nsat = 0
  character(len=32), allocatable :: nome_sat(:)
  real(dp), allocatable :: gms(:)      ! G*m_i em R_p^3/dia^2
end module sistema

!-----------------------------------------------------------------------
module forcas
  use constantes
  use sistema
  implicit none
contains

  ! Aceleracao por unidade de massa devida a J2 e J4 do planeta
  ! (raio de referencia = 1 nas unidades internas).
  subroutine acel_achatamento(r, a)
    real(dp), intent(in)  :: r(3)
    real(dp), intent(out) :: a(3)
    real(dp) :: r2, rr, mu2, f2, f4

    r2  = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
    rr  = sqrt(r2)
    mu2 = r(3)*r(3)/r2

    f2 = -1.5_dp*j2*gmp/(r2*r2*rr)                 ! -3/2 J2 GM / r^5
    f4 = 0.625_dp*j4*gmp/(r2*r2*r2*rr)             !  5/8 J4 GM / r^7

    a(1) = r(1)*( f2*(1.0_dp - 5.0_dp*mu2) + f4*(3.0_dp - 42.0_dp*mu2 + 63.0_dp*mu2*mu2) )
    a(2) = r(2)*( f2*(1.0_dp - 5.0_dp*mu2) + f4*(3.0_dp - 42.0_dp*mu2 + 63.0_dp*mu2*mu2) )
    a(3) = r(3)*( f2*(3.0_dp - 5.0_dp*mu2) + f4*(15.0_dp - 70.0_dp*mu2 + 63.0_dp*mu2*mu2) )
  end subroutine acel_achatamento

  ! Equacoes de movimento planetocentricas (inclui termos indiretos e a
  ! reacao do planeta ao achatamento). x = posicoes (3*nsat).
  ! Chamada pelo RA15 com a mesma assinatura do FORCE original.
  subroutine force(tm, x, v, f)
    real(dp), intent(in)  :: tm
    real(dp), intent(in)  :: x(:), v(:)
    real(dp), intent(out) :: f(:)
    real(dp) :: ri3(nsat), aobl(3, nsat), ind(3), d(3), d2, d3
    integer  :: i, j, ki, kj

    ! termos que dependem so de cada satelite
    ind = 0.0_dp
    do i = 1, nsat
      ki = 3*(i-1)
      ri3(i) = (x(ki+1)**2 + x(ki+2)**2 + x(ki+3)**2)**1.5_dp
      call acel_achatamento(x(ki+1:ki+3), aobl(:, i))
      ! aceleracao indireta (menos a aceleracao do planeta): o planeta e'
      ! puxado pelos satelites e sente a reacao do proprio achatamento
      ind = ind - gms(i)*x(ki+1:ki+3)/ri3(i) + (gms(i)/gmp)*aobl(:, i)
    end do

    do i = 1, nsat
      ki = 3*(i-1)
      ! direto (planeta + achatamento) + indireto (inclui o termo j = i,
      ! que da' o -(G M + G m_i) r/r^3 do problema de dois corpos)
      f(ki+1:ki+3) = -gmp*x(ki+1:ki+3)/ri3(i) + aobl(:, i) + ind
      ! perturbacao direta dos outros satelites
      do j = 1, nsat
        if (j == i .or. gms(j) == 0.0_dp) cycle
        kj = 3*(j-1)
        d  = x(kj+1:kj+3) - x(ki+1:ki+3)
        d2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
        d3 = d2*sqrt(d2)
        f(ki+1:ki+3) = f(ki+1:ki+3) + gms(j)*d/d3
      end do
    end do
  end subroutine force

  ! Energia total (vezes G) no referencial baricentrico. So' serve para
  ! monitorar a qualidade da integracao.
  function energia(x, v) result(en)
    real(dp), intent(in) :: x(:), v(:)
    real(dp) :: en, vp(3), vb(3), mtot, r, z2, p2, p4, d(3)
    integer  :: i, j, ki, kj

    mtot = gmp + sum(gms)
    vp = 0.0_dp
    do i = 1, nsat
      ki = 3*(i-1)
      vp = vp - gms(i)*v(ki+1:ki+3)
    end do
    vp = vp/mtot

    en = 0.5_dp*gmp*sum(vp*vp)
    do i = 1, nsat
      ki = 3*(i-1)
      vb = v(ki+1:ki+3) + vp
      en = en + 0.5_dp*gms(i)*sum(vb*vb)
      r  = sqrt(sum(x(ki+1:ki+3)**2))
      z2 = (x(ki+3)/r)**2
      p2 = 0.5_dp*(3.0_dp*z2 - 1.0_dp)
      p4 = (35.0_dp*z2*z2 - 30.0_dp*z2 + 3.0_dp)/8.0_dp
      en = en - gmp*gms(i)/r*(1.0_dp - j2*p2/r**2 - j4*p4/r**4)
      do j = i+1, nsat
        kj = 3*(j-1)
        d  = x(kj+1:kj+3) - x(ki+1:ki+3)
        en = en - gms(i)*gms(j)/sqrt(sum(d*d))
      end do
    end do
  end function energia

end module forcas

!-----------------------------------------------------------------------
module elementos
  use constantes
  implicit none
contains

  ! Solucao da equacao de Kepler (eliptica)
  function kepler(am, e) result(u)
    real(dp), intent(in) :: am, e
    real(dp) :: u, du, m
    integer  :: it
    m = modulo(am, pi2)
    u = m + e*sin(m)
    do it = 1, 100
      du = (u - e*sin(u) - m)/(1.0_dp - e*cos(u))
      u  = u - du
      if (abs(du) < 1.0e-15_dp) exit
    end do
  end function kepler

  ! Elementos osculadores -> estado cartesiano (angulos em rad)
  subroutine osc_para_cart(mu, a, e, inc, w, om, am, r, v)
    real(dp), intent(in)  :: mu, a, e, inc, w, om, am
    real(dp), intent(out) :: r(3), v(3)
    real(dp) :: u, cu, su, xo, yo, vxo, vyo, n, fac
    real(dp) :: p(3), q(3), cw, sw, co, so, ci, si

    u  = kepler(am, e)
    cu = cos(u);  su = sin(u)
    n  = sqrt(mu/a**3)
    fac = sqrt(1.0_dp - e*e)

    xo  = a*(cu - e)
    yo  = a*fac*su
    vxo = -n*a*su/(1.0_dp - e*cu)
    vyo =  n*a*fac*cu/(1.0_dp - e*cu)

    cw = cos(w);   sw = sin(w)
    co = cos(om);  so = sin(om)
    ci = cos(inc); si = sin(inc)

    p = [ co*cw - so*sw*ci,  so*cw + co*sw*ci, sw*si ]
    q = [-co*sw - so*cw*ci, -so*sw + co*cw*ci, cw*si ]

    r = xo*p  + yo*q
    v = vxo*p + vyo*q
  end subroutine osc_para_cart

  ! Estado cartesiano -> elementos osculadores (angulos em rad, [0,2pi))
  ! ierr = 1 se a orbita nao for eliptica (e >= 1)
  subroutine cart_para_osc(mu, r, v, a, e, inc, w, om, am, ierr)
    real(dp), intent(in)  :: mu, r(3), v(3)
    real(dp), intent(out) :: a, e, inc, w, om, am
    integer,  intent(out) :: ierr
    real(dp) :: h(3), hh, rr, v2, rv, esinu, ecosu, u, f, cosu_lat, sinu_lat, ulat

    ierr = 0
    rr = sqrt(sum(r*r))
    v2 = sum(v*v)
    rv = sum(r*v)
    h  = [ r(2)*v(3) - r(3)*v(2), r(3)*v(1) - r(1)*v(3), r(1)*v(2) - r(2)*v(1) ]
    hh = sqrt(sum(h*h))

    a = 1.0_dp/(2.0_dp/rr - v2/mu)
    if (a <= 0.0_dp) then
      ierr = 1
      e = 1.0_dp + sqrt(max(0.0_dp, 1.0_dp - hh*hh/(mu*abs(a))))
      inc = acos(max(-1.0_dp, min(1.0_dp, h(3)/hh)))
      w = 0.0_dp; om = 0.0_dp; am = 0.0_dp
      return
    end if

    esinu = rv/sqrt(mu*a)
    ecosu = rr*v2/mu - 1.0_dp
    e = sqrt(esinu*esinu + ecosu*ecosu)
    u = atan2(esinu, ecosu)
    am = modulo(u - e*sin(u), pi2)
    f  = 2.0_dp*atan2(sqrt(1.0_dp + e)*sin(0.5_dp*u), sqrt(1.0_dp - e)*cos(0.5_dp*u))

    inc = acos(max(-1.0_dp, min(1.0_dp, h(3)/hh)))

    if (sqrt(h(1)**2 + h(2)**2) < 1.0e-14_dp*hh) then
      ! orbita equatorial: nodo indefinido, usa longitude verdadeira
      om   = 0.0_dp
      ulat = atan2(r(2), r(1))
    else
      om = modulo(atan2(h(1), -h(2)), pi2)
      cosu_lat = r(1)*cos(om) + r(2)*sin(om)
      sinu_lat = (r(2)*cos(om) - r(1)*sin(om))*cos(inc) + r(3)*sin(inc)
      ulat = atan2(sinu_lat, cosu_lat)
    end if
    w = modulo(ulat - f, pi2)
  end subroutine cart_para_osc

  !---------------------------------------------------------------------
  ! Elementos geometricos (Renner & Sicardy 2006, Cel. Mech. 94, 237)
  ! Adaptado da subrotina "geo" original. Entrada: estado x(6)
  ! planetocentrico; gm, rp, cj2, cj4, cj6 nas mesmas unidades de x.
  ! Saida: a, e, dincl, w, om, am (angulos em rad, [0,2pi)).
  !---------------------------------------------------------------------
  subroutine geo(gm, rp, cj2, cj4, cj6, x, a, e, dincl, w, om, am)
    real(dp), intent(in)  :: gm, rp, cj2, cj4, cj6, x(6)
    real(dp), intent(out) :: a, e, dincl, w, om, am
    integer,  parameter   :: itmax = 1000
    real(dp), parameter   :: tol = 1.0e-13_dp        ! tolerancia relativa em a
    real(dp) :: r, L, rpt, Lpt, rc, Lc, z, zc, rptc, Lptc, zptc, zpt
    real(dp) :: n, ni, a0, inc, kapa, kapa2, eta2, qui2, c1, c2
    real(dp) :: alfa1, alfa2, alfaqu, adn, adkapa, adni, adeta2, adqui2
    real(dp) :: dt1, dt2, dt3, dt4, dt5, dt6, dt7, dt8, dt10, dt9, dt11
    real(dp) :: lambda, wtil
    integer  :: it

    z   = x(3)
    zpt = x(6)
    r   = sqrt(x(1)*x(1) + x(2)*x(2))
    L   = modulo(atan2(x(2), x(1)), pi2)

    rpt = x(4)*cos(L) + x(5)*sin(L)
    Lpt = (-x(4)*sin(L) + x(5)*cos(L))/r

    a = r;  e = 0.0_dp;  inc = 0.0_dp
    rc = 0.0_dp; Lc = 0.0_dp; zc = 0.0_dp
    rptc = 0.0_dp; Lptc = 0.0_dp; zptc = 0.0_dp
    a0 = a

    do it = 1, itmax
      c1 = sqrt(gm/a**3)
      c2 = (rp/a)**2

      n = c1*(1.0_dp + 0.75_dp*c2*cj2 - 9.0_dp/32.0_dp*c2**2*cj2**2 + &
              27.0_dp/128.0_dp*c2**3*cj2**3 + 3.0_dp*c2*cj2*e*e - 12.0_dp*c2*cj2*inc*inc)
      adn = -15.0_dp*c2*c2*cj4/16.0_dp + 45.0_dp*c2**3*cj2*cj4/64.0_dp + 35.0_dp*c2**3*cj6/32.0_dp
      n = n + c1*adn

      kapa = c1*(1.0_dp - 0.75_dp*c2*cj2 - 9.0_dp/32.0_dp*c2*c2*cj2*cj2 - &
                 27.0_dp/128.0_dp*c2**3*cj2**3 - 9.0_dp*c2*cj2*inc*inc)
      adkapa = 45.0_dp*c2*c2*cj4/16.0_dp + 135.0_dp*c2**3*cj2*cj4/64.0_dp - 175.0_dp*c2**3*cj6/32.0_dp
      kapa = kapa + c1*adkapa

      ni = c1*(1.0_dp + 2.25_dp*c2*cj2 - 81.0_dp/32.0_dp*c2*c2*cj2*cj2 + &
               729.0_dp/128.0_dp*c2**3*cj2**3 + 6.0_dp*c2*cj2*e*e - 51.0_dp/4.0_dp*c2*cj2*inc*inc)
      adni = -75.0_dp*c2*c2*cj4/16.0_dp + 675.0_dp*c2**3*cj2*cj4/64.0_dp + 245.0_dp*c2**3*cj6/32.0_dp
      ni = ni + c1*adni

      eta2   = c1*c1*(1.0_dp - 2.0_dp*c2*cj2)
      adeta2 = 75.0_dp*c2*c2*cj4/8.0_dp - 175.0_dp*c2**3*cj6/8.0_dp
      eta2   = eta2 + c1*c1*adeta2

      qui2   = c1*c1*(1.0_dp + 7.5_dp*c2*cj2)
      adqui2 = -175.0_dp*c2*c2*cj4/8.0_dp + 735.0_dp*c2**3*cj6/16.0_dp
      qui2   = qui2 + c1*c1*adqui2

      alfa1  = (2.0_dp*ni + kapa)/3.0_dp
      alfa2  = 2.0_dp*ni - kapa
      alfaqu = alfa1*alfa2
      kapa2  = kapa*kapa

      dt1 = r - rc
      dt2 = (Lpt - Lptc - n)/(n + n)
      a   = dt1/(1.0_dp - dt2)

      dt3 = (rpt - rptc)/(a*kapa)
      e   = sqrt(dt2**2 + dt3**2)

      dt4 = (z - zc)/a
      dt5 = (zpt - zptc)/(a*ni)
      inc = sqrt(dt4**2 + dt5**2)

      lambda = L - Lc - 2.0_dp*n*dt3/kapa

      dt6 = a*kapa*(1.0_dp - (r - rc)/a)
      dt7 = rpt - rptc
      dt8 = modulo(atan2(dt7, dt6), pi2)
      wtil = lambda - dt8

      dt10 = ni*(z - zc)
      dt9  = zpt - zptc
      dt11 = modulo(atan2(dt10, dt9), pi2)
      om   = lambda - dt11

      lambda = modulo(lambda, pi2)
      wtil   = modulo(wtil, pi2)
      om     = modulo(om, pi2)
      w      = modulo(wtil - om, pi2)
      am     = modulo(lambda - wtil, pi2)
      dincl  = inc

      if (abs(a - a0) <= tol*abs(a)) return
      a0 = a

      rc = a*e*e*(1.5_dp*eta2/kapa2 - 1.0_dp - 0.5_dp*eta2*cos(2.0_dp*lambda - 2.0_dp*wtil)/kapa2) + &
           a*inc*inc*(0.75_dp*qui2/kapa2 - 1.0_dp + 0.25_dp*qui2*cos(2.0_dp*lambda - 2.0_dp*om)/alfaqu)

      Lc = e*e*(0.75_dp + 0.5_dp*eta2/kapa2)*n/kapa*sin(2.0_dp*lambda - 2.0_dp*wtil) - &
           inc*inc*0.25_dp*qui2*n*sin(2.0_dp*lambda - 2.0_dp*om)/(alfaqu*ni)

      zc = a*inc*e*(0.5_dp*qui2*sin(2.0_dp*lambda - wtil - om)/(kapa*alfa1) - &
                    1.5_dp*qui2*sin(wtil - om)/(kapa*alfa2))

      rptc = a*e*e*eta2*sin(2.0_dp*lambda - 2.0_dp*wtil)/kapa - &
             a*inc*inc*0.5_dp*qui2*ni*sin(2.0_dp*lambda - 2.0_dp*om)/alfaqu

      Lptc = e*e*n*(3.5_dp - 3.0_dp*eta2/kapa2 - 0.5_dp*kapa2/(n*n) + &
             (1.5_dp + eta2/kapa2)*cos(2.0_dp*lambda - 2.0_dp*wtil)) + &
             inc*inc*n*(2.0_dp - 0.5_dp*kapa2/(n*n) - 1.5_dp*qui2/kapa2 - &
             0.5_dp*qui2/alfaqu*cos(2.0_dp*lambda - 2.0_dp*om))

      zptc = a*inc*e*(0.5_dp*qui2*(kapa + ni)*cos(2.0_dp*lambda - wtil - om)/(kapa*alfa1) + &
             1.5_dp*qui2*(kapa - ni)*cos(wtil - om)/(kapa*alfa2))
    end do
  end subroutine geo

end module elementos

!-----------------------------------------------------------------------
module integrador
  ! RA15 - Everhart (1985), mesma rotina do codigo original, apenas
  ! convertida para formato livre e com vetores de tamanho variavel.
  ! ierr = 1 se nao conseguiu escolher o passo inicial.
  use constantes
  use forcas, only: force
  implicit none
contains

  subroutine ra15(tm, x, v, tf, ll, xl, nv, nclass, ierr)
    integer,  intent(in)    :: ll, nv, nclass
    real(dp), intent(inout) :: tm, x(nv), v(nv), xl
    real(dp), intent(in)    :: tf
    integer,  intent(out)   :: ierr
    real(dp) :: f1(nv), fj(nv), y(nv), z(nv)
    real(dp) :: b(7,nv), g(7,nv), e(7,nv), bd(7,nv)
    real(dp) :: c(21), d(21), r(21), w(7), u(7)
    integer, parameter :: nw(8) = [0,0,1,3,6,10,15,21]
    real(dp), parameter :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp, sr = 1.4_dp
    real(dp), parameter :: h(8) = [0.0_dp, .05626256053692215_dp, .18024069173689236_dp, &
         .35262471711316964_dp, .54715362633055538_dp, .73421017721541053_dp, &
         .88532094683909577_dp, .97752061356128750_dp]
    logical  :: nsf, nper, npq, ncl, nes
    real(dp) :: dir, pw, ww, w1, ss, tp, t, t2, tval, s, q, a, temp, gk, hv, tval1
    integer  :: n, k, l, la, lb, lc, ld, le, ncount, ns, nf, ni, m, j, jd

    ierr = 0
    nper = .false.
    nsf  = .false.
    ncl  = nclass == 1
    npq  = nclass < 2
    nes  = ll < 0

    dir = merge(-one, one, tf < zero)
    xl  = dir*abs(xl)
    pw  = 1.0_dp/9.0_dp

    do n = 2, 8
      ww = n + n*n
      if (ncl) ww = n
      w(n-1) = one/ww
      ww = n
      u(n-1) = one/ww
    end do

    do k = 1, nv
      if (ncl) v(k) = zero
      bd(:,k) = zero
      b(:,k)  = zero
    end do

    w1 = half
    if (ncl) w1 = one
    c(1) = -h(2)
    d(1) =  h(2)
    r(1) = one/(h(3) - h(2))
    la = 1
    lc = 1
    do k = 3, 7
      lb = la
      la = lc + 1
      lc = nw(k+1)
      c(la) = -h(k)*c(lb)
      c(lc) = c(la-1) - h(k)
      d(la) = h(2)*d(lb)
      d(lc) = -c(lc)
      r(la) = one/(h(k+1) - h(2))
      r(lc) = one/(h(k+1) - h(k))
      if (k /= 3) then
        do l = 4, k
          ld = la + l - 3
          le = lb + l - 4
          c(ld) = c(le) - h(k)*c(le+1)
          d(ld) = d(le) + h(l-1)*d(le+1)
          r(ld) = one/(h(k+1) - h(l-1))
        end do
      end if
    end do

    ss = 10.0_dp**(-ll)

    if (nes) then
      tp = xl
    else if (xl /= zero) then
      tp = xl
    else
      tp = 0.1_dp*dir
    end if
    if (tp/tf > half) tp = half*tf
    ncount = 0

    ! ---- inicio (ou reinicio com passo menor)
1000 continue
    ns = 0
    nf = 0
    ni = 6
    tm = zero
    call force(tm, x, v, f1)
    nf = nf + 1

    ! ---- laco de passos
2000 continue
    do k = 1, nv
      g(1,k) = b(1,k) + d(1)*b(2,k) + d(2)*b(3,k) + d(4)*b(4,k) + d(7)*b(5,k) + d(11)*b(6,k) + d(16)*b(7,k)
      g(2,k) =          b(2,k) + d(3)*b(3,k) + d(5)*b(4,k) + d(8)*b(5,k) + d(12)*b(6,k) + d(17)*b(7,k)
      g(3,k) =                   b(3,k) + d(6)*b(4,k) + d(9)*b(5,k) + d(13)*b(6,k) + d(18)*b(7,k)
      g(4,k) =                            b(4,k) + d(10)*b(5,k) + d(14)*b(6,k) + d(19)*b(7,k)
      g(5,k) =                                     b(5,k) + d(15)*b(6,k) + d(20)*b(7,k)
      g(6,k) =                                              b(6,k) + d(21)*b(7,k)
      g(7,k) =                                                       b(7,k)
    end do

    t  = tp
    t2 = t*t
    if (ncl) t2 = t
    tval = abs(t)

    do m = 1, ni
      do j = 2, 8
        jd = j - 1
        s  = h(j)
        q  = s
        if (ncl) q = one
        do k = 1, nv
          a = w(3)*b(3,k) + s*(w(4)*b(4,k) + s*(w(5)*b(5,k) + s*(w(6)*b(6,k) + s*w(7)*b(7,k))))
          y(k) = x(k) + q*(t*v(k) + t2*s*(f1(k)*w1 + s*(w(1)*b(1,k) + s*(w(2)*b(2,k) + s*a))))
          if (.not. npq) then
            a = u(3)*b(3,k) + s*(u(4)*b(4,k) + s*(u(5)*b(5,k) + s*(u(6)*b(6,k) + s*u(7)*b(7,k))))
            z(k) = v(k) + s*t*(f1(k) + s*(u(1)*b(1,k) + s*(u(2)*b(2,k) + s*a)))
          end if
        end do

        call force(tm + s*t, y, z, fj)
        nf = nf + 1

        do k = 1, nv
          temp = g(jd,k)
          gk = (fj(k) - f1(k))/s
          select case (j)
          case (2)
            g(1,k) = gk
          case (3)
            g(2,k) = (gk - g(1,k))*r(1)
          case (4)
            g(3,k) = ((gk - g(1,k))*r(2) - g(2,k))*r(3)
          case (5)
            g(4,k) = (((gk - g(1,k))*r(4) - g(2,k))*r(5) - g(3,k))*r(6)
          case (6)
            g(5,k) = ((((gk - g(1,k))*r(7) - g(2,k))*r(8) - g(3,k))*r(9) - g(4,k))*r(10)
          case (7)
            g(6,k) = (((((gk - g(1,k))*r(11) - g(2,k))*r(12) - g(3,k))*r(13) - g(4,k))*r(14) &
                     - g(5,k))*r(15)
          case (8)
            g(7,k) = ((((((gk - g(1,k))*r(16) - g(2,k))*r(17) - g(3,k))*r(18) - g(4,k))*r(19) &
                     - g(5,k))*r(20) - g(6,k))*r(21)
          end select

          temp = g(jd,k) - temp
          b(jd,k) = b(jd,k) + temp

          select case (j)
          case (3)
            b(1,k) = b(1,k) + c(1)*temp
          case (4)
            b(1,k) = b(1,k) + c(2)*temp
            b(2,k) = b(2,k) + c(3)*temp
          case (5)
            b(1,k) = b(1,k) + c(4)*temp
            b(2,k) = b(2,k) + c(5)*temp
            b(3,k) = b(3,k) + c(6)*temp
          case (6)
            b(1,k) = b(1,k) + c(7)*temp
            b(2,k) = b(2,k) + c(8)*temp
            b(3,k) = b(3,k) + c(9)*temp
            b(4,k) = b(4,k) + c(10)*temp
          case (7)
            b(1,k) = b(1,k) + c(11)*temp
            b(2,k) = b(2,k) + c(12)*temp
            b(3,k) = b(3,k) + c(13)*temp
            b(4,k) = b(4,k) + c(14)*temp
            b(5,k) = b(5,k) + c(15)*temp
          case (8)
            b(1,k) = b(1,k) + c(16)*temp
            b(2,k) = b(2,k) + c(17)*temp
            b(3,k) = b(3,k) + c(18)*temp
            b(4,k) = b(4,k) + c(19)*temp
            b(5,k) = b(5,k) + c(20)*temp
            b(6,k) = b(6,k) + c(21)*temp
          end select
        end do
      end do

      if (.not. nes .or. m >= ni) then
        hv = zero
        do k = 1, nv
          hv = max(hv, abs(b(7,k)))
        end do
        tval1 = tval**7
        hv = hv*w(7)/tval1
      end if
    end do

    if (.not. nsf) then
      if (nes) then
        tp = xl
      else
        tp = (ss**pw)/(hv**pw)*dir
        if (tp/t <= one) then
          tp = 0.8_dp*tp
          ncount = ncount + 1
          if (ncount > 10) then
            ierr = 1
            return
          end if
          goto 1000
        end if
      end if
      nsf = .true.
    end if

    do k = 1, nv
      x(k) = x(k) + v(k)*t + t2*(f1(k)*w1 + b(1,k)*w(1) + b(2,k)*w(2) + b(3,k)*w(3) + &
             b(4,k)*w(4) + b(5,k)*w(5) + b(6,k)*w(6) + b(7,k)*w(7))
      if (.not. ncl) then
        v(k) = v(k) + t*(f1(k) + b(1,k)*u(1) + b(2,k)*u(2) + b(3,k)*u(3) + &
               b(4,k)*u(4) + b(5,k)*u(5) + b(6,k)*u(6) + b(7,k)*u(7))
      end if
    end do

    tm = tm + t
    ns = ns + 1

    if (nper) return
    call force(tm, x, v, f1)
    nf = nf + 1

    if (nes) then
      tp = xl
    else
      tp = dir*(ss**pw)/(hv**pw)
      if (tp/t > sr) tp = t*sr
    end if

    if (dir*(tm + tp) >= dir*tf - 1.0e-8_dp) then
      tp = tf - tm
      nper = .true.
    end if

    q = tp/t
    do k = 1, nv
      if (ns /= 1) then
        do j = 1, 7
          bd(j,k) = b(j,k) - e(j,k)
        end do
      end if
      e(1,k) = q*(b(1,k) + 2.0_dp*b(2,k) + 3.0_dp*b(3,k) + 4.0_dp*b(4,k) + 5.0_dp*b(5,k) + &
               6.0_dp*b(6,k) + 7.0_dp*b(7,k))
      e(2,k) = q**2*(b(2,k) + 3.0_dp*b(3,k) + 6.0_dp*b(4,k) + 10.0_dp*b(5,k) + 15.0_dp*b(6,k) + &
               21.0_dp*b(7,k))
      e(3,k) = q**3*(b(3,k) + 4.0_dp*b(4,k) + 10.0_dp*b(5,k) + 20.0_dp*b(6,k) + 35.0_dp*b(7,k))
      e(4,k) = q**4*(b(4,k) + 5.0_dp*b(5,k) + 15.0_dp*b(6,k) + 35.0_dp*b(7,k))
      e(5,k) = q**5*(b(5,k) + 6.0_dp*b(6,k) + 21.0_dp*b(7,k))
      e(6,k) = q**6*(b(6,k) + 7.0_dp*b(7,k))
      e(7,k) = q**7*b(7,k)
      do l = 1, 7
        b(l,k) = e(l,k) + bd(l,k)
      end do
    end do

    ni = 2
    goto 2000
  end subroutine ra15

end module integrador

!-----------------------------------------------------------------------
module leitura
  ! Leitura do arquivo .in e do banco de dados de planetas
  use constantes
  use sistema
  implicit none

  ! parametros da simulacao
  character(len=256) :: arq_planetas = 'planetas.dat'
  character(len=256) :: pasta_saida  = '.'
  character(len=16)  :: unidade_massa = 'kg'
  real(dp) :: tempo_final = -1.0_dp, passo_saida = -1.0_dp
  integer  :: precisao = 12
  logical  :: saida_osc = .true., saida_geo = .true.
  ! elementos iniciais lidos (km, graus) e massas como no arquivo
  real(dp), allocatable :: massa_in(:), el_in(:,:)
  ! valores do planeta dados diretamente no .in (sobrescrevem o banco)
  real(dp) :: gm_in = -1.0_dp, raio_in = -1.0_dp, j2_in = -999.0_dp, j4_in = -999.0_dp

contains

  function minusculo(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i, c
    t = s
    do i = 1, len(s)
      c = iachar(s(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) t(i:i) = achar(c + 32)
    end do
  end function minusculo

  ! remove comentario (depois de '!' ou '#') e espacos
  function limpa(linha) result(s)
    character(len=*), intent(in) :: linha
    character(len=len(linha)) :: s
    integer :: k
    s = linha
    k = scan(s, '!#')
    if (k > 0) s(k:) = ' '
    do k = 1, len(s)
      if (s(k:k) == achar(9)) s(k:k) = ' '     ! tab -> espaco
    end do
    s = adjustl(s)
  end function limpa

  logical function sim_ou_nao(valor)
    character(len=*), intent(in) :: valor
    select case (minusculo(trim(adjustl(valor))))
    case ('sim', 's', 'yes', 'y', '1', '.true.', 'true', 't')
      sim_ou_nao = .true.
    case default
      sim_ou_nao = .false.
    end select
  end function sim_ou_nao

  subroutine erro(msg)
    character(len=*), intent(in) :: msg
    write(*,'(a)') ' ERRO: '//trim(msg)
    stop 1
  end subroutine erro

  subroutine le_entrada(arquivo)
    character(len=*), intent(in) :: arquivo
    character(len=512) :: linha, s, chave, valor
    integer :: u, ios, k, nmax
    character(len=32) :: nome
    real(dp) :: el(7)
    character(len=32), allocatable :: tmp_nome(:)
    real(dp), allocatable :: tmp_el(:,:)

    open(newunit=u, file=trim(arquivo), status='old', action='read', iostat=ios)
    if (ios /= 0) call erro('nao consegui abrir o arquivo de entrada "'//trim(arquivo)//'"')

    nmax = 100
    allocate(tmp_nome(nmax), tmp_el(7, nmax))
    nsat = 0

    do
      read(u, '(a)', iostat=ios) linha
      if (ios /= 0) exit
      s = limpa(linha)
      if (len_trim(s) == 0) cycle

      k = index(s, '=')
      if (k > 0) then
        chave = minusculo(trim(adjustl(s(1:k-1))))
        valor = trim(adjustl(s(k+1:)))
        select case (trim(chave))
        case ('planeta')
          nome_planeta = valor
        case ('arquivo_planetas')
          arq_planetas = valor
        case ('tempo_final')
          read(valor, *, iostat=ios) tempo_final
        case ('passo_saida')
          read(valor, *, iostat=ios) passo_saida
        case ('precisao')
          read(valor, *, iostat=ios) precisao
        case ('saida_osculadores')
          saida_osc = sim_ou_nao(valor)
        case ('saida_geometricos')
          saida_geo = sim_ou_nao(valor)
        case ('pasta_saida')
          pasta_saida = valor
        case ('unidade_massa')
          unidade_massa = minusculo(valor)
        case ('gm_planeta')
          read(valor, *, iostat=ios) gm_in
        case ('raio_planeta')
          read(valor, *, iostat=ios) raio_in
        case ('j2')
          read(valor, *, iostat=ios) j2_in
        case ('j4')
          read(valor, *, iostat=ios) j4_in
        case default
          call erro('chave desconhecida no .in: "'//trim(chave)//'"')
        end select
        if (ios /= 0) call erro('valor invalido para "'//trim(chave)//'": '//trim(valor))
        cycle
      end if

      if (minusculo(trim(s)) == 'satelites:' .or. minusculo(trim(s)) == 'satelites') then
        ! bloco de satelites ate' "fim"
        do
          read(u, '(a)', iostat=ios) linha
          if (ios /= 0) exit
          s = limpa(linha)
          if (len_trim(s) == 0) cycle
          if (minusculo(trim(s)) == 'fim') exit
          read(s, *, iostat=ios) nome, el
          if (ios /= 0) call erro('linha de satelite invalida: '//trim(s))
          nsat = nsat + 1
          if (nsat > nmax) call erro('mais de 100 satelites; aumente nmax em le_entrada')
          tmp_nome(nsat) = nome
          tmp_el(:, nsat) = el
        end do
        cycle
      end if

      call erro('linha nao reconhecida no .in: '//trim(s))
    end do
    close(u)

    if (nsat == 0) call erro('nenhum satelite no bloco "satelites:" ... "fim"')
    if (tempo_final <= 0.0_dp) call erro('defina tempo_final (em dias)')
    if (passo_saida <= 0.0_dp) call erro('defina passo_saida (em dias)')

    allocate(nome_sat(nsat), massa_in(nsat), el_in(6, nsat), gms(nsat))
    nome_sat = tmp_nome(1:nsat)
    massa_in = tmp_el(1, 1:nsat)
    el_in    = tmp_el(2:7, 1:nsat)
  end subroutine le_entrada

  ! Procura o planeta em arq_planetas. Formato de cada linha:
  !   nome[/apelido/...]   GM(km^3/s^2)   R_eq(km)   J2   J4
  subroutine le_planeta()
    character(len=512) :: linha, s
    character(len=128) :: nomes, alvo, um
    real(dp) :: vals(4)
    integer  :: u, ios, k
    logical  :: achou

    achou = .false.
    alvo = minusculo(trim(adjustl(nome_planeta)))

    if (len_trim(alvo) > 0 .and. alvo /= 'personalizado') then
      open(newunit=u, file=trim(arq_planetas), status='old', action='read', iostat=ios)
      if (ios /= 0) call erro('nao consegui abrir "'//trim(arq_planetas)//'"')
      do
        read(u, '(a)', iostat=ios) linha
        if (ios /= 0) exit
        s = limpa(linha)
        if (len_trim(s) == 0) cycle
        ! (nao da' para ler o nome com read(*): '/' encerra a leitura)
        k = index(s, ' ')
        nomes = s(1:k-1)
        read(s(k:), *, iostat=ios) vals
        if (ios /= 0) cycle
        nomes = minusculo(nomes)
        ! compara com cada apelido separado por '/'
        do
          k = index(nomes, '/')
          if (k > 0) then
            um = nomes(1:k-1)
            nomes = nomes(k+1:)
          else
            um = nomes
          end if
          if (trim(um) == trim(alvo)) then
            achou = .true.
            exit
          end if
          if (k == 0) exit
        end do
        if (achou) then
          gm_planeta_km = vals(1)
          raio_km = vals(2)
          j2 = vals(3)
          j4 = vals(4)
          exit
        end if
      end do
      close(u)
      if (.not. achou) call erro('planeta "'//trim(nome_planeta)//'" nao encontrado em '//trim(arq_planetas))
    end if

    ! valores dados no .in tem prioridade
    if (gm_in   > 0.0_dp)    gm_planeta_km = gm_in
    if (raio_in > 0.0_dp)    raio_km = raio_in
    if (j2_in   > -998.0_dp) j2 = j2_in
    if (j4_in   > -998.0_dp) j4 = j4_in

    if (gm_planeta_km <= 0.0_dp .or. raio_km <= 0.0_dp) &
      call erro('planeta sem GM/raio: use "planeta = Nome" ou gm_planeta/raio_planeta/j2/j4')
    if (len_trim(nome_planeta) == 0) nome_planeta = 'personalizado'
  end subroutine le_planeta

end module leitura

!=======================================================================
program nbody
  use constantes
  use sistema
  use forcas
  use elementos
  use integrador
  use leitura
  implicit none

  character(len=256) :: arq_in, arq
  real(dp), allocatable :: x(:), v(:), rr(:,:), vv(:,:)
  integer,  allocatable :: u_osc(:), u_geo(:)
  logical,  allocatable :: avisado(:)
  real(dp) :: conv_gm, t, tm, xl, e0, e1, emax, en
  real(dp) :: a, e, inc, w, om, am, st(6), mu
  integer  :: i, k, nv, npassos, ierr, iel

  if (command_argument_count() >= 1) then
    call get_command_argument(1, arq_in)
  else
    arq_in = 'simulacao.in'
  end if

  call le_entrada(arq_in)
  call le_planeta()

  ! --- conversao de unidades: km^3/s^2 -> R_p^3/dia^2
  conv_gm = dia*dia/raio_km**3
  gmp = gm_planeta_km*conv_gm

  do i = 1, nsat
    select case (trim(unidade_massa))
    case ('kg')
      gms(i) = G_SI_km*massa_in(i)*conv_gm
    case ('razao')                      ! m / M_planeta
      gms(i) = massa_in(i)*gmp
    case ('gm')                         ! G*m em km^3/s^2
      gms(i) = massa_in(i)*conv_gm
    case default
      call erro('unidade_massa deve ser kg, razao ou gm')
    end select
  end do

  ! --- condicoes iniciais (elementos osculadores planetocentricos)
  nv = 3*nsat
  allocate(x(nv), v(nv), rr(3,nsat), vv(3,nsat), avisado(nsat))
  avisado = .false.
  do i = 1, nsat
    k = 3*(i-1)
    mu = gmp + gms(i)
    call osc_para_cart(mu, el_in(1,i)/raio_km, el_in(2,i), el_in(3,i)*deg2rad, &
                       el_in(4,i)*deg2rad, el_in(5,i)*deg2rad, el_in(6,i)*deg2rad, &
                       x(k+1:k+3), v(k+1:k+3))
  end do

  ! --- resumo na tela
  write(*,'(a)')        ' ================= NBODY ================='
  write(*,'(a,a)')      ' Arquivo de entrada : ', trim(arq_in)
  write(*,'(a,a)')      ' Planeta            : ', trim(nome_planeta)
  write(*,'(a,es20.12)')' GM (km^3/s^2)      : ', gm_planeta_km
  write(*,'(a,f14.4)')  ' Raio equat. (km)   : ', raio_km
  write(*,'(a,es20.12)')' J2                 : ', j2
  write(*,'(a,es20.12)')' J4                 : ', j4
  write(*,'(a,i0)')     ' Satelites          : ', nsat
  do i = 1, nsat
    write(*,'(3x,a16,a,es12.5,a,f14.3,a)') nome_sat(i), '  m/M =', gms(i)/gmp, &
      '   a =', el_in(1,i), ' km'
  end do
  write(*,'(a,f16.4,a)') ' Tempo final        : ', tempo_final, ' dias'
  write(*,'(a,f16.8,a)') ' Passo de saida     : ', passo_saida, ' dias'
  write(*,'(a,i0)')      ' Precisao (LL RA15) : ', precisao

  ! --- abre arquivos de saida
  allocate(u_osc(nsat), u_geo(nsat))
  do i = 1, nsat
    if (saida_osc) then
      arq = trim(pasta_saida)//'/'//trim(nome_sat(i))//'_osc.dat'
      open(newunit=u_osc(i), file=trim(arq), status='replace')
      write(u_osc(i),'(a)') '# Elementos OSCULADORES planetocentricos - '//trim(nome_sat(i))// &
                            ' (planeta: '//trim(nome_planeta)//')'
      write(u_osc(i),'(a)') '#  tempo(dias)          a(km)                 e                  ' // &
                            'inc(deg)       w(deg)        om(deg)       M(deg)'
    end if
    if (saida_geo) then
      arq = trim(pasta_saida)//'/'//trim(nome_sat(i))//'_geo.dat'
      open(newunit=u_geo(i), file=trim(arq), status='replace')
      write(u_geo(i),'(a)') '# Elementos GEOMETRICOS (Renner & Sicardy 2006) - '//trim(nome_sat(i))// &
                            ' (planeta: '//trim(nome_planeta)//')'
      write(u_geo(i),'(a)') '#  tempo(dias)          a(km)                 e                  ' // &
                            'inc(deg)       w(deg)        om(deg)       M(deg)'
    end if
  end do

  ! --- integracao
  e0 = energia(x, v)
  emax = 0.0_dp
  t = 0.0_dp
  call escreve_saida()

  npassos = nint(tempo_final/passo_saida)
  do iel = 1, npassos
    xl = 0.0_dp
    ! NCLASS = -2: y'' = F(y), sem dependencia na velocidade
    call ra15(tm, x, v, passo_saida, precisao, xl, nv, -2, ierr)
    if (ierr /= 0) then
      write(*,'(a,f16.6)') ' ERRO: RA15 nao conseguiu escolher o passo em t = ', t
      exit
    end if
    t = iel*passo_saida

    ! colisao com o planeta?
    do i = 1, nsat
      k = 3*(i-1)
      if (sum(x(k+1:k+3)**2) < 1.0_dp) then
        write(*,'(a,a,a,f16.6)') ' AVISO: ', trim(nome_sat(i)), ' colidiu com o planeta em t = ', t
        call fecha_e_para()
      end if
    end do

    call escreve_saida()
    if (e0 /= 0.0_dp) then
      en = energia(x, v)
      emax = max(emax, abs((en - e0)/e0))
    end if

    if (mod(iel, max(1, npassos/20)) == 0) &
      write(*,'(a,f7.1,a,f16.4,a)') '   ', 100.0_dp*iel/npassos, '%   t =', t, ' dias'
  end do

  if (e0 /= 0.0_dp) then
    e1 = energia(x, v)
    write(*,'(a,es10.3)') ' Erro relativo final da energia : ', abs((e1 - e0)/e0)
    write(*,'(a,es10.3)') ' Erro relativo maximo da energia: ', emax
  end if
  call fecha_e_para()

contains

  subroutine escreve_saida()
    integer :: j, kk, ie
    do j = 1, nsat
      kk = 3*(j-1)
      st(1:3) = x(kk+1:kk+3)
      st(4:6) = v(kk+1:kk+3)
      if (saida_osc) then
        call cart_para_osc(gmp + gms(j), st(1:3), st(4:6), a, e, inc, w, om, am, ie)
        if (ie /= 0 .and. .not. avisado(j)) then
          write(*,'(a,a,a,f16.6)') ' AVISO: orbita de ', trim(nome_sat(j)), &
                                   ' nao e mais eliptica em t = ', t
          avisado(j) = .true.
        end if
        write(u_osc(j), 100) t, a*raio_km, e, inc*rad2deg, w*rad2deg, om*rad2deg, am*rad2deg
      end if
      if (saida_geo) then
        call geo(gmp, 1.0_dp, j2, j4, 0.0_dp, st, a, e, inc, w, om, am)
        write(u_geo(j), 100) t, a*raio_km, e, inc*rad2deg, w*rad2deg, om*rad2deg, am*rad2deg
      end if
    end do
100 format(f18.6, 1x, es21.13, 1x, es21.13, 4(1x, f13.8))
  end subroutine escreve_saida

  subroutine fecha_e_para()
    integer :: j
    do j = 1, nsat
      if (saida_osc) close(u_osc(j))
      if (saida_geo) close(u_geo(j))
    end do
    write(*,'(a)') ' Fim.'
    stop
  end subroutine fecha_e_para

end program nbody
