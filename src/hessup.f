

C      hessup.fpp ---> hessup.f
C     update the hessian via BFGS
C
C     compile as:
C-R      g77 -E -P -DUSING_R -o hessup.f hessup.fpp
C     or
C-S+     g77 -E -P  -o hessup.f hessup.fpp
C     prior to make or INSTALL

      subroutine hessup( dgr, dparm, newgr, oldgr, curprm, oldprm,
     $     reset, bk, newhss, bks, qrbk, paradj, np ,newprm)

      integer i,j,reset,np,iq1,iqy1,iq100,rank,piv(30),dx(2)

      double precision dgr(np), dparm(np), newgr(np),
     $     oldgr(np), curprm(np), oldprm(np), bk(np,np),
     $     newhss(np,np), bks(np), qrbk(np,np), paradj(np),
     $     tmp, dgrdpr, dprbks, tol,wrksp(60),newprm(np)

      data tol /1.0D-09/

C###			d.grad <- as.matrix(new$grad - old$grad)
C#			d.parm <- as.matrix(cur.parm - old.parm)

  
      do i = 1,np
         dgr(i) = newgr(i) - oldgr(i)
         dparm(i) = curprm(i) - oldprm(i)
      end do
      
C#	if(reset.hess) {
      
      if (reset.eq.1) then
C#		bk <- new$hess
         do i = 1, np
            do j = 1,np
            bk(i,j) = newhss(i,j)
         end do
      end do
C     #		reset.hess <- FALSE
      reset = 0
      

C#	}
      end if 
      
      

C#			bks <- bk %*% d.parm
C     also get denoms for use in later loop

      dgrdpr = 0.0
      dprbks = 0.0

      
      do i = 1,np
         tmp=0.0
         do j = 1,np
            tmp = tmp + bk(i,j)*dparm(j)
         end do
         bks(i) = tmp
         dgrdpr = dgrdpr + dgr(i)*dparm(i)
         dprbks = dprbks + dparm(i)*bks(i)
      end do


C#			bk <- bk - bks %*% t(bks)/c(t(d.parm) %*% bks) + d.grad %*%
C#				t(d.grad)/c(t(d.grad) %*% d.parm)

      do i = 1,np
         do j = 1,np
            tmp = dgr(i)*dgr(j)/dgrdpr - bks(i)*bks(j)/dprbks 
            bk(i,j) = bk(i,j) + tmp
         end do
      end do
C#            cur.parm <- cur.parm - qr.coef(qr(bk), new$grad)

      
      do i = 1,np
         do j = 1,np
            qrbk(i,j)=bk(i,j)
         end do
      end do
      do i = 1, np
         piv(i)=i
      end do

      
      rank = np
      dx(1)=np
      dx(2)=np




      call dqrdc2(qrbk,np,np,np,tol,rank,newhss,piv,wrksp)

      iq1=1
      iqy1=1
      iq100=100

          




      call dqrsl(qrbk,np,np,np,newhss,newgr,wrksp,paradj,paradj,wrksp,
     $     wrksp,100,iq1)


      do i = 1,np
         curprm(i) = curprm(i) - paradj(i)
      end do


      return
      end

