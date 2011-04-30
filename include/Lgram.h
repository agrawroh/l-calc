/*

   Copyright (C) Michael Rubinstein

   This file is part of the L-function package L.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   Check the License for details. You should have received a copy of it, along
   with the package; see the file 'COPYING'. If not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

*/


    //returns first gram point >t
    template <class ttype>
    Double L_function <ttype>:: initialize_gram(Double t)
    {

        Double N1=Nmain(t)/2,N2;
        Double t2=t;
        Long n1,n2;
        Double incr=1.01/(Nmain(t+1)/2-N1);
        Double x,tmp;
        Double y;
        Long u;

        if (incr<0)incr=.1;

        n1=(Long)(floor(lcalc_to_double(N1))+.1);

        //Nmain(T) behaves funny near the origin. So, if we are
        //asking for a small gram point, start at 20, and scan backwards
        if(N1<=2.){
            t2=12;
            do{
                t2=t2-.2;
                tmp=Nmain(t2)/2;
            }while(t2>t&&tmp>0);
            if(t2<t) {
                t2=t;
                //if zeta, dirichlet, cusp form for SL2(Z), or maass L-function, 
                //we can return 0 as a gram point
                if(t==0&&what_type_L<4&&what_type_L!=0) return 0;
            }
            else {n1=-1;incr=.2;}



        }


            //fudge with the increment until we're sure to capture the next gram point
            N2=Nmain(t2+incr)/2;
            n2=(Long)(floor(lcalc_to_double(N2))+.1); if(N2<0)n2--;
            do{
                if((n2-n1)<1){incr=incr*1.4; N2=Nmain(t2+incr)/2;n2=(Long)(floor(lcalc_to_double(N2))+.1);if(N2<0)n2--;}
                if((n2-n1)>=2){incr=incr*.9; N2=Nmain(t2+incr)/2;n2=(Long)(floor(lcalc_to_double(N2))+.1);if(N2<0)n2--;}
               //cout << t2 << " " << t2+incr << " " << n1 << " " << n2 << endl;
            }while((n2-n1)<1||(n2-n1)>=2);

            y=t2+incr;
            //divide and conquer to compute the next gram point
            for(int i=1;i<=20;i++){
                 x=(t2+y)/2;  tmp=Nmain(x)/2;
                 u=(Long) (floor(lcalc_to_double(tmp))+.1); if(tmp<0) u--;
                 if(u==n1){ t2=x;}
                 else{ y=x;}
                 //cout << t2 <<  " " << y << endl;
            }

        return x;

    }

    //computes the gram point above the current gram point t, i.e.
    //assumes that t is itself a gram point.
    //I don't want to use initialize_gram(current_gram) because
    //the numeric current gram might be slightly smaller than actual current
    //gram and we might just get the same point
    template <class ttype>
    Double L_function <ttype>:: next_gram(Double t)
    {
        Double N1=(Nmain(t))/2,N2;
        Double t2=t;
        Long n1,n2;
        Double incr=1.01/((Nmain(t+1))/2-N1);

        if (incr<0)incr=.1;

        n1=(Long)(rint(lcalc_to_double(N1))+.1); //n1 is the closest integer to Nmain(t)/2


        //fudge with the increment until we're sure to capture the next gram point
        N2=(Nmain(t2+incr))/2;
        n2=(Long)(floor(lcalc_to_double(N2))+.1);
        do{
            if((n2-n1)<1){incr=incr*1.4; N2=Nmain(t2+incr)/2;n2=(Long)(floor(lcalc_to_double((N2)))+.1);}
            if((n2-n1)>=2){incr=incr*.9; N2=Nmain(t2+incr)/2;n2=(Long)(floor(lcalc_to_double((N2)))+.1);}
           //cout << N1 << " " << N2 << " " << n1 << " " << n2 << endl;
        }while((n2-n1)<1||(n2-n1)>=2);

        Double x;
        Double y=t2+incr;
        //divide and conquer to compute the next gram point
        for(int i=1;i<=20;i++){
             x=(t2+y)/2;
             if(((Long)(floor(lcalc_to_double(Nmain(x)/2))+.1))==n1){ t2=x;}
             else{ y=x;}
             //cout << t2 <<  " " << y << endl;
        }
        //cout << x << " " << Nmain(x)/2 << endl;
        //cout << "========================\n";
        return x;
    }

    //find the nth gram point t such that Nmain(t)/2=n
    template <class ttype>
    Double L_function <ttype>:: nth_gram(Long n)
    {


         if(n==0) return 0.;

         Double t=1000.,r,x,b=Double(n)*Double(n)*tolerance_sqrd;
         do{
             x=(n-.5*Nmain(t));
             r=.5*(Nmain(t-.05)-Nmain(t+.05))*10; //this is essentially Newton's method, but with a crude approximation
                                        //for the derivative. Works reasonably well. Reason for approx: in case
                                        //n, hence t, is large, we won't have much precision left over for the derivative.
             t-=x/r;
             if(my_verbose>3)
                 cout << "nth_gram("<< n << "): N(" << t << " )/2=" << Nmain(t)/2 << "  , difference:" << x << endl;
         }while(x*x>b);

         return t;
    }

