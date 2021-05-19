function Jmnp = LH_INTEGRALS_OLD(a,r,z,alpha,K,E,cas)
%LH_INTEGRALS Lipschitz-Hankel integrals.
%   Jmnp = LH_INTEGRALS(r,z,alpha,K,E,cas) calculates the Jmnp LH integral
%   at the point (r,z) and at a cut angle 'alpha' (see Paynter et al., 
%   2009).
%
%   The integral is characterized by the case 'cas', passed as a string
%   of the integral index 'mnp'. For example, J001 at
%   a point (0,1) is:
%       LH_INTEGRALS(0,1,0,K,E,'001')
%
%   'K' and 'E' are the complete elliptic integrals of the first and
%   kind for the given point (r,z). It is calculated outside of the
%   function to improve computational time.
%
%   University of Oxford 
%   Department of Engineering Science
%   Jhonatan Da Ponte Lopes, MSc 
%   May, 2017; Last revision: 2017-05-24

%-----------------------------------------------------------------------
%                            VERSION NUMBER
%-----------------------------------------------------------------------

%-------|--------------|----------------------------------|------------|
% V.N.  | Modified by  | Modifications made               | Date       |
%-------|--------------|----------------------------------|------------|
% 1.00  | J.P. Lopes   | File creation                    | 2017-05-24 |
% 1.10  | J.P. Lopes   | Add LH Integrals - 02P           | 2018-11-08 |


%-------------------------------------------------------------------
%                              CHECKS
%-------------------------------------------------------------------

par=a;

%[~,~,P]=elliptic123(m,n)

% Parameters
k2=(4.*par.*r)./((par+r).^2+z.^2);
kp2=1-k2;
h=(4.*par.*r)./((par+r).^2);
%P = EllipticPi(h,sqrt(k2));

switch cas

%-------------------------------------------------------------------
%                              J_00P
%-------------------------------------------------------------------

    case '000'
        if r==0
            Jmnp=sqrt(1./(par.^2+z.^2));
        else
            Jmnp=1/pi.*sqrt(k2./(par.*r)).*K;
        end

    case '001' %OK
        Jmnp=(2.*z.*E)./(pi.*((par-r).^2+z.^2).*sqrt((par+r).^2+z.^2));
        
    case '002' %OK
        Jmnp=-((2.*(((par.^2-r.^2).^2-2.*(par.^2+r.^2).*z.^2-3.*z.^4).*...
            E+z.^2.*((par-r).^2+z.^2).*K))./(pi.*((par-r).^2+z.^2).^2.*...
            ((par+r).^2+z.^2).^(3./2)));
        
    case '003' %OK
        Jmnp=(2.*z.*((-12.*(par.^2-r.^2).^2.*(par.^2+r.^2)+...
            (-13.*par.^4+58.*par.^2.*r.^2-13.*r.^4).*z.^2+...
            10.*(par.^2+r.^2).*z.^4+11.*z.^6).*E+((par-r).^2+z.^2).*...
            (3.*(par.^2-r.^2).^2-2.*(par.^2+r.^2).*z.^2-5.*z.^4).*K))./...
            (pi.*((par-r).^2+z.^2).^3.*((par+r).^2+z.^2).^(5./2));

%-------------------------------------------------------------------
%                              J_01P
%-------------------------------------------------------------------

    case '01m1'
        if r>par
            corr2=1;
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %h=(4.*par.*r)./((par+r).^2);
            P = EllipticPi(h,sqrt(k2));
            %P = elliptic3(pi/2,k2,h);
            %P = ellipticPi(k2,h);
            corr1=(r-par)./(r+par).*P;
        elseif r==par
            corr2=1/2;
            corr1=0;
        else
            corr2=0;
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %h=(4.*par.*r)./((par+r).^2);
            P = EllipticPi(h,sqrt(k2));
            %P = elliptic3(pi/2,k2,h);
            %P = ellipticPi(k2,h);
            corr1=((-par+r).*P)./(par+r);
        end
        Jmnp=2/pi.*(sqrt(par./(r.*k2)).*E)+(r.^2-par.^2)./...
            (pi.*r.*sqrt((par.^2+r.^2).^2+z.^2)).*K+...
            z.^2./(2.*pi.*r).*sqrt(k2./(r.*par)).*corr1-...
            z./r.*sign(-alpha+atan(z./(-par+r))).*corr2;

    case '010' %OK
        if r>par
            corr2=1./r;
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %h=(4.*par.*r)./((par+r).^2);
            P = EllipticPi(h,sqrt(k2));
            %P = elliptic3(pi/2,k2,h);
            %P = ellipticPi(k2,h);
            corr1=((-par+r).*P)./(par+r);
        elseif r==par
            corr2=1./(2.*par);
            corr1=0;
        elseif r==0
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %h=(4.*par.*r)./((par+r).^2);
            P = EllipticPi(h,sqrt(k2));
            %P = elliptic3(pi/2,k2,h);
            %P = ellipticPi(k2,h);
            corr1=(-par+r)./(par+r).*P;
            corr2=0;
        else
            corr2=0;
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %h=(4.*par.*r)./((par+r).^2);
            P = EllipticPi(h,sqrt(k2));
            %P = elliptic3(pi/2,k2,h);
            %P = ellipticPi(k2,h);
            corr1=((-par+r).*P)./(par+r);
        end
        Jmnp=(z.*sqrt(1./((par+r).^2+z.^2)).*(-K-corr1))./...
            (pi.*r)+corr2.*sign(-alpha+atan(z./(-par+r)));
        
    case '010byr' %OK
        if r==0
            Jmnp=z./2.*(1./(par.^2+z.^2)).^(3/2);
        else
            if r>par
                corr2=1./r;
                %k2=(4.*par.*r)./((par+r).^2+z.^2);
                %h=(4.*par.*r)./((par+r).^2);
                P = EllipticPi(h,sqrt(k2));
                %P = elliptic3(pi/2,k2,h);
                %P = ellipticPi(k2,h);
                corr1=((-par+r).*P)./(par+r);
            elseif r==par
                corr2=1./(2.*par);
                corr1=0;
            else
                corr2=0;
                %k2=(4.*par.*r)./((par+r).^2+z.^2);
                %h=(4.*par.*r)./((par+r).^2);
                P = EllipticPi(h,sqrt(k2));
                %P = elliptic3(pi/2,k2,h);
                %P = ellipticPi(k2,h);
                corr1=((-par+r).*P)./(par+r);
            end
            Jmnp=1./(r).*((z.*sqrt(1./((par+r).^2+z.^2)).*(-K-corr1))./...
                (pi.*r)+corr2.*sign(-alpha+atan(z./(-par+r))));
        end
        
    case '011' %OK
        if r==0
            Jmnp=0;
        else
            Jmnp=1./(2.*pi.*r).*sqrt(k2./(par.*r)).*(K-(par.^2-r.^2+z.^2)./...
                ((par-r).^2+z.^2).*E);
        end
        
    case '011byr' %OK
        if r==0
            Jmnp=-(1./2).*(par.^2-2.*z.^2).*(1./(par.^2+z.^2)).^(5/2);
        else
            Jmnp=(1./(2.*pi.*r).*sqrt(k2./(par.*r)).*(K-(par.^2-r.^2+z.^2)./...
                ((par-r).^2+z.^2).*E))./r;
        end
        
    case '012' %OK
        if r==0
            Jmnp=0;
        else
            Jmnp=-((z.*((par.^4-7.*r.^4-6.*r.^2.*z.^2+z.^4+2.*par.^2.*...
                (3.*r.^2+z.^2)).*E-((par-r).^2+z.^2).*(par.^2-r.^2+z.^2).*...
                K))./(pi.*r.*((par-r).^2+z.^2).^2.*((par+r).^2+z.^2).^(3./2)));
        end
        
    case '012byr' %OK
        if r==0
            Jmnp=3/2.*z.*(1/(par.^2+z.^2)).^(7/2).*(-3.*par.^2+2.*z.^2);
        else
            Jmnp=-((z.*((par.^4-7.*r.^4-6.*r.^2.*z.^2+z.^4+2.*par.^2.*...
                (3.*r.^2+z.^2)).*E-((par-r).^2+z.^2).*(par.^2-r.^2+z.^2).*...
                K))./(pi.*r.*((par-r).^2+z.^2).^2.*((par+r).^2+z.^2).^(3./2)))./(r);
        end
        
    case '013' %OK
        if r==0
            Jmnp=0;
        else
            Jmnp=(((par.^2-r.^2).^3.*(par.^2+7.*r.^2)+...
                (par.^6-75.*par.^4.*r.^2+51.*par.^2.*r.^4+23.*r.^6).*z.^2+...
                (-3.*par.^4-46.*par.^2.*r.^2+65.*r.^4).*z.^4+...
                (-5.*par.^2+33.*r.^2).*z.^6-2.*z.^8).*E-((par-r).^2+z.^2).*...
                ((par.^2-r.^2).^3+12.*r.^2.*(-par+r).*(par+r).*z.^2+...
                (-3.*par.^2+11.*r.^2).*z.^4-2.*z.^6).*K)./...
                (pi.*r.*((par-r).^2+z.^2).^3.*((par+r).^2+z.^2).^(5./2));
        end

%-------------------------------------------------------------------
%                              J_02P
%-------------------------------------------------------------------

    case '020'
        if r>par
            corr2=1./r;
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %h=(4.*par.*r)./((par+r).^2);
            P = EllipticPi(h,sqrt(k2));
            %P = elliptic3(pi/2,k2,h);
            %P = ellipticPi(k2,h);
            corr1=((-par+r).*P)./(par+r);
        elseif r==par
            corr2=1./(2.*par);
            corr1=0;
        else
            corr2=0;
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %h=(4.*par.*r)./((par+r).^2);
            P = EllipticPi(h,sqrt(k2));
            %P = elliptic3(pi/2,k2,h);
            %P = ellipticPi(k2,h);
            corr1=((-par+r).*P)./(par+r);
        end
        Jmnp=-1./pi.*sqrt(k2./(par.*r)).*(-4.*par./(r.*k2).*E+...
            par.^2./(r.^2).*K-z.^2./(r.^2).*corr1)-...
            2.*z./(r).*corr2.*sign(-alpha+atan(z./(-par+r)));
        
        
    case '021'
        if r>par
            corr2=1./r;
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %h=(4.*par.*r)./((par+r).^2);
            P = EllipticPi(h,sqrt(k2));
            %P = elliptic3(pi/2,k2,h);
            %P = ellipticPi(k2,h);
            corr1=((-par+r).*P)./(par+r);
        elseif r==par
            corr2=1./(2.*par);
            corr1=0;
        else
            corr2=0;
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %h=(4.*par.*r)./((par+r).^2);
            P = EllipticPi(h,sqrt(k2));
            %P = elliptic3(pi/2,k2,h);
            %P = ellipticPi(k2,h);
            corr1=((-par+r).*P)./(par+r);
        end
        Jmnp=-2.*z./(pi.*r.^2.*((par.^2-r.^2).^2+z.^2).*...
            sqrt((par.^2+r.^2).^2+z.^2)).*(r.^2.*E+...
            ((par.^2-r.^2).^2+z.^2).*(K+corr1))+...
            2./r.*corr2.*sign(-alpha+atan(z./(-par+r)));
       

%-------------------------------------------------------------------
%                              J_10P
%-------------------------------------------------------------------
    
    case '10m1' %OK
        if r<par
            corr2=1;
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %h=(4.*par.*r)./((par+r).^2);
            P = EllipticPi(h,sqrt(k2));
            %P = elliptic3(pi/2,k2,h);
            %P = ellipticPi(k2,h);
            corr1=((par-r).*P)./(par+r);
        elseif r==par
            corr2=1/2;
            corr1=0;
        else
            corr2=0;
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %h=(4.*par.*r)./((par+r).^2);
            P = EllipticPi(h,sqrt(k2));
            %P = elliptic3(pi/2,k2,h);
            %P = ellipticPi(k2,h);
            corr1=((par-r).*P)./(par+r);
        end
        Jmnp=(sqrt(((par+r).^2+z.^2)./(par.^2)).*E)./pi+...
            ((par.^2-r.^2).*K)./(pi.*par.*sqrt((par+r).^2+z.^2))+...
            (z.^2.*sqrt(1./((par+r).^2+z.^2)).*corr1)./(pi.*par)-...
            (z.*corr2.*sign(-alpha+atan(z./(par-r))))./(par);

    case '100' %OK
        if r<par
            corr2=1./par;
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %h=(4.*par.*r)./((par+r).^2);
            P = EllipticPi(h,sqrt(k2));
            %P = elliptic3(pi/2,k2,h);
            %P = ellipticPi(k2,h);
            corr1=((par-r).*P)./(par+r);
        elseif r==par
            corr2=1./(2.*r);
            corr1=0;
        else
            corr2=0;
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %h=(4.*par.*r)./((par+r).^2);
            P = EllipticPi(h,sqrt(k2));
            %P = elliptic3(pi/2,k2,h);
            %P = ellipticPi(k2,h);
            corr1=((par-r).*P)./(par+r);
        end
        Jmnp=(z.*sqrt(1./((par+r).^2+z.^2)).*(-K-corr1))./...
            (pi.*par)+corr2.*sign(-alpha+atan(z./(par-r)));
        

    case '101' %OK
        Jmnp=(sqrt(1./((par+r).^2+z.^2)).*(-(((-par.^2+r.^2+z.^2).*E)./...
            ((-par+r).^2+z.^2))+K))./(par.*pi);
        
    case '102' %OK
        Jmnp=-((z.*((-7.*par.^4+6.*par.^2.*(r-z).*(r+z)+...
            (r.^2+z.^2).^2).*E+(par.^2-r.^2-z.^2).*...
            ((par-r).^2+z.^2).*K))./(par.*pi.*((par-r).^2+z.^2).^2.*...
            ((par+r).^2+z.^2).^(3/2)));
    
    case '103' %OK
        Jmnp=((-7.*par.^8+(r.^2-2.*z.^2).*(r.^2+z.^2).^3+...
            par.^6.*(20.*r.^2+23.*z.^2)+...
            par.^2.*(r.^2+z.^2).*(4.*r.^4-79.*r.^2.*z.^2+33.*z.^4)+...
            par.^4.*(-18.*r.^4+51.*r.^2.*z.^2+65.*z.^4)).*E+...
            ((par-r).^2+z.^2).*((par.^2-r.^2).^3+12.*par.^2.*...
            (-par.^2+r.^2).*z.^2+(-11.*par.^2+3.*r.^2).*z.^4+...
            2.*z.^6).*K)./(par.*pi.*((par-r).^2+z.^2).^3.*...
            ((par+r).^2+z.^2).^(5/2));
        
%-------------------------------------------------------------------
%                              J_11P
%-------------------------------------------------------------------

    case '11m1' %OK
        if r==0
            Jmnp=0;
        else
            J010=LH_INTEGRALS(a,r,z,alpha,K,E,'010');
            J100=LH_INTEGRALS(a,r,z,alpha,K,E,'100');
            J110=LH_INTEGRALS(a,r,z,alpha,K,E,'110');
            Jmnp=1./2.*(par.*J010+r.*J100-z.*J110);
        end
        
    case '11m1byr' %OK
        if r==0
            Jmnp=1/2.*(-(z./sqrt(1+z.^2))+sign(-alpha+atan(z)));
        else
            J010byr=LH_INTEGRALS(a,r,z,alpha,K,E,'010byr');
            J100=LH_INTEGRALS(a,r,z,alpha,K,E,'100');
            J110byr=LH_INTEGRALS(a,r,z,alpha,K,E,'110byr');
            Jmnp=1./2.*(par.*J010byr+J100-z.*J110byr);
        end
        
    case '110' %OK
        if r==0
            Jmnp=0;
        else
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            Jmnp=2./(pi.*sqrt(par.*r.*k2)).*((2-k2)./2.*K-E);
        end

    case '110byr' %OK
        if r==0
            Jmnp=((par.^2./(par.^2+z.^2)).^(3/2))./(2.*par.^2);
        else
            Jmnp=(-E+1/2.*(2-(4.*par.*r)./((par+r).^2+z.^2)).*K)./...
                (pi.*r.*sqrt((par.^2.*r.^2)./((par+r).^2+z.^2)));
        end
        
    case '111' %OK
        if r==0
            Jmnp=0;
        else
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %kp2=1-k2;
            Jmnp=z./(par.*pi.*r.*sqrt((par+r).^2+z.^2)).*...
                ((2-k2)./(2.*kp2).*E-K);
        end
        
    case '111byr' %OK
        if r==0
            Jmnp=(3.*par.*z)./(2.*(par.^2+z.^2).^(5/2));
        else
            %k2=(4.*par.*r)./((par+r).^2+z.^2);
            %kp2=1-k2;
            Jmnp=1./(r).*(z./(par.*pi.*r.*sqrt((par+r).^2+z.^2)).*...
                ((2-k2)./(2.*kp2).*E-K));
        end
        
    case '112' %OK
        if r==0
            Jmnp=0;
        else
            Jmnp=(-(par.^6-par.^4.*(r.^2-2.*z.^2)+r.^2.*(r.^2+z.^2).^2+...
                par.^2.*(-r.^4-12.*r.^2.*z.^2+z.^4)).*E+(par.^6-...
                2.*par.^5.*r-par.^4.*(r.^2-2.*z.^2)-2.*par.*r.^3.*...
                (r.^2+z.^2)+r.^2.*(r.^2+z.^2).^2+par.^3.*...
                (4.*r.^3-2.*r.*z.^2)+par.^2.*(-r.^4+z.^4)).*K)./...
                (par.*pi.*r.*(par.^2-2.*par.*r+r.^2+z.^2).^2.*...
                ((par+r).^2+z.^2).^(3/2));
        end
        
    case '113' %OK
        if r==0
            Jmnp=0;
        else
            Jmnp=-((z.*((3.*par.^8+3.*r.^2.*(r.^2+z.^2).^3+9.*par.^6.*...
                (4.*r.^2+z.^2)+par.^4.*(-78.*r.^4-41.*r.^2.*z.^2+9.*z.^4)+...
                par.^2.*(36.*r.^6-41.*r.^4.*z.^2-74.*r.^2.*z.^4+3.*z.^6)).*E+...
                (-3.*par.^8+6.*par.^7.*r-9.*par.^6.*z.^2+6.*par.*r.^3.*...
                (r.^2+z.^2).^2-3.*r.^2.*(r.^2+z.^2).^3-6.*par.^5.*...
                (r.^3-2.*r.*z.^2)+par.^4.*(6.*r.^4+17.*r.^2.*z.^2-9.*z.^4)+...
                par.^3.*(-6.*r.^5-40.*r.^3.*z.^2+6.*r.*z.^4)+...
                par.^2.*(17.*r.^4.*z.^2+14.*r.^2.*z.^4-3.*z.^6)).*K))./...
                (par.*pi.*r.*(par.^2-2.*par.*r+r.^2+z.^2).^3.*...
                ((par+r).^2+z.^2).^(5/2)));
        end

%-------------------------------------------------------------------
%                              J_12P
%-------------------------------------------------------------------        

    case '12m1'
        % Not checked for r==0!
        J020=LH_INTEGRALS(a,r,z,alpha,K,E,'020');
        J110=LH_INTEGRALS(a,r,z,alpha,K,E,'100');
        J120=LH_INTEGRALS(a,r,z,alpha,K,E,'100');
        Jmnp=1/3.*(par.*J020+r.*J110-z.*J120);
        
    case '120'
        % Not checked for r==0!
        J11m1=LH_INTEGRALS(a,r,z,alpha,K,E,'11m1');
        J100=LH_INTEGRALS(a,r,z,alpha,K,E,'100');
        Jmnp=2./r.*J11m1-J100;
        
    case '121'
        % Not checked for r==0!
        J101=LH_INTEGRALS(a,r,z,alpha,K,E,'101');
        J110=LH_INTEGRALS(a,r,z,alpha,K,E,'110');
        Jmnp=2./r.*J110-J101;
        
    case '122'
        % Not checked for r==0!
        J103=LH_INTEGRALS(a,r,z,alpha,K,E,'101');
        J112=LH_INTEGRALS(a,r,z,alpha,K,E,'110');
        Jmnp=2./r.*J112-J103;

%-------------------------------------------------------------------
%                              J_20P
%-------------------------------------------------------------------

    case '20m1' %OK
        J11m1=LH_INTEGRALS(a,r,z,alpha,K,E,'11m1');
        J10m1=LH_INTEGRALS(a,r,z,alpha,K,E,'10m1');
        Jmnp=-r./(par).*J11m1-z./(par).*J10m1;

    case '201' %OK
        J100=LH_INTEGRALS(a,r,z,alpha,K,E,'100');
        J001=LH_INTEGRALS(a,r,z,alpha,K,E,'001');
        Jmnp=2./(par).*J100-J001;
        
    case '202' %OK
        J101=LH_INTEGRALS(a,r,z,alpha,K,E,'101');
        J002=LH_INTEGRALS(a,r,z,alpha,K,E,'002');
        Jmnp=2./(par).*J101-J002;
        
    case '203' %OK
        J102=LH_INTEGRALS(a,r,z,alpha,K,E,'102');
        J003=LH_INTEGRALS(a,r,z,alpha,K,E,'003');
        Jmnp=2./(par).*J102-J003;
        
        
%-------------------------------------------------------------------
%                              J_21P
%-------------------------------------------------------------------

    case '210'%(NOT DEFINED FOR r==0!) %OK
        J11m1=LH_INTEGRALS(a,r,z,alpha,K,E,'11m1');
        J010=LH_INTEGRALS(a,r,z,alpha,K,E,'010');
        Jmnp=2./(par).*J11m1-J010;
        
    case '210byr' %OK
        J11m1byr=LH_INTEGRALS(a,r,z,alpha,K,E,'11m1byr');
        J010byr=LH_INTEGRALS(a,r,z,alpha,K,E,'010byr');
        Jmnp=2./(par).*J11m1byr-J010byr;
        
    case '211' %OK
        J110=LH_INTEGRALS(a,r,z,alpha,K,E,'110');
        J011=LH_INTEGRALS(a,r,z,alpha,K,E,'011');
        Jmnp=2./(par).*J110-J011;
        
    case '211byr' %OK
        if r==0
            Jmnp=(par.^3.*sqrt(1./(par.^2+z.^2))-2.*par.*z.^2.*...
                sqrt(1./(par.^2+z.^2))+2.*par.^2.*...
                sqrt(par.^2./(par.^2+z.^2))+2.*z.^2.*...
                sqrt(par.^2./(par.^2+z.^2)))./(2.*par.*(par.^2+z.^2).^2);
        else
            J110byr=LH_INTEGRALS(a,r,z,alpha,K,E,'110byr');
            J011byr=LH_INTEGRALS(a,r,z,alpha,K,E,'011byr');
            Jmnp=2./(par).*J110byr-J011byr;
        end
        
    case '212' %OK
        J111=LH_INTEGRALS(a,r,z,alpha,K,E,'111');
        J012=LH_INTEGRALS(a,r,z,alpha,K,E,'012');
        Jmnp=2./(par).*J111-J012;
        
    case '212byr' %OK
        if r==0
            Jmnp=3/2.*z.*(-2.*(1./(par.^2+z.^2)).^(5/2)+5.*par.^2.*...
                (1./(par.^2+z.^2)).^(7/2)+2./(par.^2+z.^2).^(5/2));
        else
            J212=LH_INTEGRALS(a,r,z,alpha,K,E,'212');
            Jmnp=1./(r).*J212;
        end
        
    case '213' %OK
        J112=LH_INTEGRALS(a,r,z,alpha,K,E,'112');
        J013=LH_INTEGRALS(a,r,z,alpha,K,E,'013');
        Jmnp=2./(par).*J112-J013;
        
    case '213byr' %OK
        if r==0
            Jmnp=3/2.*(-3.*(1./(par.^2+z.^2)).^(5/2)-...
                35.*z.^4.*(1./(par.^2+z.^2)).^(9/2)-...
                2./(par.^2+z.^2).^(5/2)+10.*z.^2.*...
                (3.*(1./(par.^2+z.^2)).^(7/2)+1./(par.^2+z.^2).^(7/2)));
        else
            J213=LH_INTEGRALS(a,r,z,alpha,K,E,'213');
            Jmnp=1./(r).*J213;
        end
                
%-------------------------------------------------------------------
%                              J_22P
%-------------------------------------------------------------------
    case '220'%(NOT DEFINED FOR r==0!)
        J12m1=LH_INTEGRALS(a,r,z,alpha,K,E,'12m1');
        J020=LH_INTEGRALS(a,r,z,alpha,K,E,'010');
        Jmnp=2./(par).*J12m1-J020;
        
    case '221'
        J120=LH_INTEGRALS(a,r,z,alpha,K,E,'120');
        J021=LH_INTEGRALS(a,r,z,alpha,K,E,'021');
        Jmnp=2./(par).*J120-J021;
        
        
    otherwise
        error('LH integral not recognised!');
end


end