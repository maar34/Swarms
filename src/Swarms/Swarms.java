package Swarms;


import com.cycling74.max.*;
import com.cycling74.jitter.*;
import java.util.Random;


public class Swarms extends MaxObject
{

	public static class Methods2D {

//		static AnimationPanel2D panel = new AnimationPanel2D( 1, 1);


	    // ============================
	    public static  void updateCeldas(Celda2D[][] C, double[] x, double[] y){
	    	int L = C.length;
	    	int N = x.length;
			int[][] nC = new int[L][L]; // numero de particulas en cada Celda2D
			int n,ix,iy;
			
			for(ix = 0; ix < L; ++ix){
			    for(iy = 0; iy < L; ++iy){
				nC[ix][iy] = 0;
			    }
			}
			
			// encuentra el numero de particulas en cada Celda2D
			for(n = 0; n < N; ++n){
			    ix = (int)x[n];
			    iy = (int)y[n];
			    ++nC[ix][iy];
			}
			// define los arreglos de cada Celda2D para guardar
			// las particulas que se encuentran en dicha Celda2D
			for(ix = 0; ix < L; ++ix){
			    for(iy = 0; iy < L; ++iy){
					C[ix][iy].setNP(nC[ix][iy]);
					nC[ix][iy] = 0;
			    }
			}
			// llena los arreglos de cada Celda2D con los indices
			// de las particulas que se encuentran en dicha Celda2D
			for(n = 0; n < N; ++n){
			    ix = (int)x[n];
			    iy = (int)y[n];
			    C[ix][iy].setPar(nC[ix][iy],n);
			    ++nC[ix][iy];
	    		
			}
	    } // Fin de updateCeldas
	    
	    //  ================================


	    public static double updatePsi(double[] a){
		
		double psix = 0;
		double psiy = 0;
		double psi;	
		int N = a.length;
		double  vx, vy;
		

		for(int n = 0; n < N; ++n){
	    vx = Math.cos(a[n]);
	    vy = Math.sin(a[n]);//CALCULO DE PSI
	    psix += vx;
	    psiy += vy;
		}

		psix /= N;
		psiy /= N;
		psi = Math.sqrt((psix*psix) + (psiy*psiy));
		
		return psi;
		}




	    public static double updateVicsek(Celda2D[][] C, double[] x, double[] y, double[] a, 
	    		double radius, double eta,double u){
	    	int L = C.length;
	    	int N = x.length;
	    	double[] angle = new double[N];
	    	int lx,ly,lxv,lyv,cx,cy;
	    	int m,n,nv,ip,ipv;
	    	double dx,dy,vx,vy,dr;
	    	double psi,psiX,psiY;
		

	    	// ciclos sobre las Celda2Ds
	    	for(lx = 0; lx < L; ++lx){
	    		for(ly = 0; ly < L; ++ly){
	    			// ciclo sobre las particulas de cada Celda2D
	    			for(n = 0; n < C[lx][ly].getNP(); ++n){
	    				ip = C[lx][ly].getPar(n);
	    				vx = 0;
	    				vy = 0;
	    				// ciclo sobre las Celda2Ds vecinas
	    				for(m = 0; m < 9; ++m){
	    					lxv = C[lx][ly].getVecinoX(m);
	    					lyv = C[lx][ly].getVecinoY(m);
	    					cx = C[lx][ly].getPhantomCompX(m);
	    					cy = C[lx][ly].getPhantomCompY(m);
	    					// ciclo sobre las particulas de la Celda2D vecina
	    					for(nv = 0; nv < C[lxv][lyv].getNP(); ++nv){
	    						ipv = C[lxv][lyv].getPar(nv);
	    						dx = Math.abs(x[ip]-x[ipv]-cx);
	    						dy = Math.abs(y[ip]-y[ipv]-cy);
	    						dr = Math.sqrt((dx*dx)+(dy*dy));
	    						//if((dx <= radius) && (dy <= radius)){
	    						if(dr <= radius){
	    							vx += Math.cos(a[ipv]);
	    							vy += Math.sin(a[ipv]);
	    						}
	    					}
	    				}
	    				angle[ip] = Math.atan2(vy,vx);
	    			}
	    		}
	    	}
	    	psi = 0;
	    	psiX = 0;
	    	psiY = 0;
	    	for(n = 0; n < N; ++n){
	    		a[n] = angle[n] + ((2*Math.PI*eta)*(Math.random() - 0.5));
	    		vx = u*Math.cos(a[n]);
	    		vy = u*Math.sin(a[n]);
	    		psiX += vx;
	    		psiY += vy;
	    		x[n] += vx;
	    		y[n] += vy;
	    		if(x[n] < 0)
	    			x[n] += L;
	    		if(x[n] >= L)
	    			x[n] -= L;
	    		if(y[n] < 0)
	    			y[n] += L;
	    		if(y[n] >= L)
	    			y[n] -= L;
	    		
	    	}
	    	psi = Math.sqrt((psiX*psiX)+(psiY*psiY));
	    	psi /= (u*N);
	    	return psi;
	    }// Fin de updateVicsek

	    // =============================
	    public static double updateChate(Celda2D[][] C,double[] x,double[] y, double[] a,
	    		double radius,double eta, double u)
	    {
	    	int L = C.length;
	    	int N = x.length;
	    	double[] vx = new double[N];
	    	double[] vy = new double[N];
	    	int lx,ly,lxv,lyv,cx,cy;
	    	int m,n,nv,ip,ipv,nvecinos;
	    	double dx,dy,da,dr,psiX,psiY;
	    	double psi;
	    	
	    	// ciclos sobre las Celda2Ds
	    	for(lx = 0; lx < L; ++lx){
	    		for(ly = 0; ly < L; ++ly){
	    			// ciclo sobre las particulas de cada Celda2D
	    			for(n = 0; n < C[lx][ly].getNP(); ++n){
	    				ip = C[lx][ly].getPar(n);
	    				vx[ip] = 0;
	    				vy[ip] = 0;
	    				nvecinos = 0;
	    				// ciclo sobre las Celda2Ds vecinas
	    				for(m = 0; m < 9; ++m){
	    					lxv = C[lx][ly].getVecinoX(m);
	    					lyv = C[lx][ly].getVecinoY(m);
	    					cx = C[lx][ly].getPhantomCompX(m);
	    					cy = C[lx][ly].getPhantomCompY(m);
	    					// 	ciclo sobre las particulas de la Celda2D vecina
	    					for(nv = 0; nv < C[lxv][lyv].getNP(); ++nv){
	    						ipv = C[lxv][lyv].getPar(nv);
	    						dx = Math.abs(x[ip]-x[ipv]-cx);
	    						dy = Math.abs(y[ip]-y[ipv]-cy);
	    						dr = Math.sqrt((dx*dx)+(dy*dy));
	    						//if((dx <= radius) && (dy <= radius)){
	    						if(dr <= radius){
	    							vx[ip] += Math.cos(a[ipv]);
	    							vy[ip] += Math.sin(a[ipv]);
	    							++nvecinos;
	    						}
	    					}
	    				}
	    				vx[ip] = vx[ip]/nvecinos;
	    				vy[ip] = vy[ip]/nvecinos;
	    			}
	    		}
	    	}
	    	psiX = 0;
	    	psiY = 0;
	    	for(n = 0; n < N; ++n){
	    		da = 2*Math.PI*Math.random();
	    		dx = vx[n] + eta*Math.cos(da);
	    		dy = vy[n] + eta*Math.sin(da);
	    		
	    		a[n] = Math.atan2(dy,dx);
	    		dx = u*Math.cos(a[n]);
	    		dy = u*Math.sin(a[n]);
	    		psiX += dx;
	    		psiY += dy;
	    		
	    		x[n] += dx;
	    		y[n] += dy;
	    		if(x[n] < 0)
	    			x[n] += L;
	    		if(x[n] >= L)
	    			x[n] -= L;
	    		if(y[n] < 0)
	    			y[n] += L;
	    		if(y[n] >= L)
	    			y[n] -= L;
	    	}
	    	psi = Math.sqrt((psiX*psiX)+(psiY*psiY));
	    	psi /= (u*N);
	    	return psi;
	    }// Fin de updateChate

	    //  =============================
	    public static double updateVicsekChate(Celda2D[][] C,double[] x,double[] y, double[] a,
	    		double radius,double eta1, double eta2,double u)
	    {
	    	int L = C.length;
	    	int N = x.length;
	    	double[] vx = new double[N];
	    	double[] vy = new double[N];
	    	int lx,ly,lxv,lyv,cx,cy;
	    	int m,n,nv,ip,ipv,nvecinos;
	    	double dx,dy,da,dr,dens;
	    	
	    	// ciclos sobre las Celda2Ds
	    	dens = 0;
	    	for(lx = 0; lx < L; ++lx){
	    		for(ly = 0; ly < L; ++ly){
	    			// ciclo sobre las particulas de cada Celda2D
	    			for(n = 0; n < C[lx][ly].getNP(); ++n){
	    				ip = C[lx][ly].getPar(n);
	    				vx[ip] = 0;
	    				vy[ip] = 0;
	    				nvecinos = 0;
	    				// ciclo sobre las Celda2Ds vecinas
	    				for(m = 0; m < 9; ++m){
	    					lxv = C[lx][ly].getVecinoX(m);
	    					lyv = C[lx][ly].getVecinoY(m);
	    					cx = C[lx][ly].getPhantomCompX(m);
	    					cy = C[lx][ly].getPhantomCompY(m);
	    					// 	ciclo sobre las particulas de la Celda2D vecina
	    					for(nv = 0; nv < C[lxv][lyv].getNP(); ++nv){
	    						ipv = C[lxv][lyv].getPar(nv);
	    						dx = Math.abs(x[ip]-x[ipv]-cx);
	    						dy = Math.abs(y[ip]-y[ipv]-cy);
	    						dr = Math.sqrt((dx*dx)+(dy*dy));
	    						//if((dx < radius) && (dy < radius)){
	    						if(dr < radius){
	    							vx[ip] += Math.cos(a[ipv]);
	    							vy[ip] += Math.sin(a[ipv]);
	    							++nvecinos;
	    						}
	    					}
	    				}
	    				vx[ip] = vx[ip]/nvecinos;
	    				vy[ip] = vy[ip]/nvecinos;
	    				dens += nvecinos;
	    			}
	    		}

	    	}
	    	dens /= N;
	    	for(n = 0; n < N; ++n){
	    		da = 2*Math.PI*Math.random();
	    		dx = vx[n] + eta2*Math.cos(da);
	    		dy = vy[n] + eta2*Math.sin(da);
	    		
	    		a[n] = Math.atan2(dy,dx);
	    		a[n] += (2*Math.PI*eta1*(Math.random()-0.5));
	    		dx = u*Math.cos(a[n]);
	    		dy = u*Math.sin(a[n]);
	    		
	    		x[n] += dx;
	    		y[n] += dy;
	    		if(x[n] < 0)
	    			x[n] += L;
	    		if(x[n] >= L)
	    			x[n] -= L;
	    		if(y[n] < 0)
	    			y[n] += L;
	    		if(y[n] >= L)
	    			y[n] -= L;
	
	    	}
	    	return dens;
	    }// Fin de updateVicsekChate
	    
	    //  ================================
	    public  static double updateVicsekRandomMixing(Celda2D[][] C, double[] x, double[] y, double[] a, 
	    		double radius, double eta){
	    	int L = C.length;
	    	int N = x.length;
	    	double[] angle = new double[N];
	    	int lx,ly,lxv,lyv,cx,cy;
	    	int m,n,nv,ip,ipv;
	    	double dx,dy,vx,vy,dr;
	    	double psi,psiX,psiY;
		

	    	// ciclos sobre las Celda2Ds
	    	for(lx = 0; lx < L; ++lx){
	    		for(ly = 0; ly < L; ++ly){
	    			// ciclo sobre las particulas de cada Celda2D
	    			for(n = 0; n < C[lx][ly].getNP(); ++n){
	    				ip = C[lx][ly].getPar(n);
	    				vx = 0;
	    				vy = 0;
	    				// ciclo sobre las Celda2Ds vecinas
	    				for(m = 0; m < 9; ++m){
	    					lxv = C[lx][ly].getVecinoX(m);
	    					lyv = C[lx][ly].getVecinoY(m);
	    					cx = C[lx][ly].getPhantomCompX(m);
	    					cy = C[lx][ly].getPhantomCompY(m);
	    					// ciclo sobre las particulas de la Celda2D vecina
	    					for(nv = 0; nv < C[lxv][lyv].getNP(); ++nv){
	    						ipv = C[lxv][lyv].getPar(nv);
	    						dx = Math.abs(x[ip]-x[ipv]-cx);
	    						dy = Math.abs(y[ip]-y[ipv]-cy);
	    						dr = Math.sqrt((dx*dx)+(dy*dy));
	    						//if((dx <= radius) && (dy <= radius)){
	    						if(dr <= radius){
	    							vx += Math.cos(a[ipv]);
	    							vy += Math.sin(a[ipv]);
	    						}
	    					}
	    				}
	    				angle[ip] = Math.atan2(vy,vx);
	    			}
	    		}
	    	}
	    	psi = 0;
	    	psiX = 0;
	    	psiY = 0;
	    	for(n = 0; n < N; ++n){
	    		a[n] = angle[n] + ((2*Math.PI*eta)*(Math.random() - 0.5));
	    		vx = Math.cos(a[n]);
	    		vy = Math.sin(a[n]);
	    		psiX += vx;
	    		psiY += vy;
	    		x[n] = L*Math.random();
	    		y[n] = L*Math.random();
	    		
	    	}
	    	psi = Math.sqrt((psiX*psiX)+(psiY*psiY));
	    	psi /= (1.0*N);
	    	return psi;
	    }// Fin de updateVicsek

	    // =============================
	    public static double updateChateRandomMixing(Celda2D[][] C,double[] x,double[] y, double[] a,
	    		double radius,double eta)
	    {
	    	int L = C.length;
	    	int N = x.length;
	    	double[] vx = new double[N];
	    	double[] vy = new double[N];
	    	int lx,ly,lxv,lyv,cx,cy;
	    	int m,n,nv,ip,ipv,nvecinos;
	    	double dx,dy,da,dr,psiX,psiY;
	    	double psi;
	    	
	    	// ciclos sobre las Celda2Ds
	    	for(lx = 0; lx < L; ++lx){
	    		for(ly = 0; ly < L; ++ly){
	    			// ciclo sobre las particulas de cada Celda2D
	    			for(n = 0; n < C[lx][ly].getNP(); ++n){
	    				ip = C[lx][ly].getPar(n);
	    				vx[ip] = 0;
	    				vy[ip] = 0;
	    				nvecinos = 0;
	    				// ciclo sobre las Celda2Ds vecinas
	    				for(m = 0; m < 9; ++m){
	    					lxv = C[lx][ly].getVecinoX(m);
	    					lyv = C[lx][ly].getVecinoY(m);
	    					cx = C[lx][ly].getPhantomCompX(m);
	    					cy = C[lx][ly].getPhantomCompY(m);
	    					// 	ciclo sobre las particulas de la Celda2D vecina
	    					for(nv = 0; nv < C[lxv][lyv].getNP(); ++nv){
	    						ipv = C[lxv][lyv].getPar(nv);
	    						dx = Math.abs(x[ip]-x[ipv]-cx);
	    						dy = Math.abs(y[ip]-y[ipv]-cy);
	    						dr = Math.sqrt((dx*dx)+(dy*dy));
	    						//if((dx <= radius) && (dy <= radius)){
	    						if(dr <= radius){
	    							vx[ip] += Math.cos(a[ipv]);
	    							vy[ip] += Math.sin(a[ipv]);
	    							++nvecinos;
	    						}
	    					}
	    				}
	    				vx[ip] = vx[ip]/nvecinos;
	    				vy[ip] = vy[ip]/nvecinos;
	    			}
	    		}
	    	}
	    	psiX = 0;
	    	psiY = 0;
	    	for(n = 0; n < N; ++n){
	    		da = 2*Math.PI*Math.random();
	    		dx = vx[n] + eta*Math.cos(da);
	    		dy = vy[n] + eta*Math.sin(da);
	    		
	    		a[n] = Math.atan2(dy,dx);
	    		dx = Math.cos(a[n]);
	    		dy = Math.sin(a[n]);
	    		psiX += dx;
	    		psiY += dy;
	    		
	    		x[n] = L*Math.random();
	    		y[n] = L*Math.random();
	    	}
	    	psi = Math.sqrt((psiX*psiX)+(psiY*psiY));
	    	psi /= (1.0*N);
	    	return psi;
	    }// Fin de updateChate

	    //  =============================
	    public static double updateVicsekChateRandomMixing(Celda2D[][] C,double[] x,double[] y, double[] a,
	    		double radius,double eta1, double eta2)
	    {
	    	int L = C.length;
	    	int N = x.length;
	    	double[] vx = new double[N];
	    	double[] vy = new double[N];
	    	int lx,ly,lxv,lyv,cx,cy;
	    	int m,n,nv,ip,ipv,nvecinos;
	    	double dx,dy,da,dr,dens;
	    	
	    	// ciclos sobre las Celda2Ds
	    	dens = 0;
	    	for(lx = 0; lx < L; ++lx){
	    		for(ly = 0; ly < L; ++ly){
	    			// ciclo sobre las particulas de cada Celda2D
	    			for(n = 0; n < C[lx][ly].getNP(); ++n){
	    				ip = C[lx][ly].getPar(n);
	    				vx[ip] = 0;
	    				vy[ip] = 0;
	    				nvecinos = 0;
	    				// ciclo sobre las Celda2Ds vecinas
	    				for(m = 0; m < 9; ++m){
	    					lxv = C[lx][ly].getVecinoX(m);
	    					lyv = C[lx][ly].getVecinoY(m);
	    					cx = C[lx][ly].getPhantomCompX(m);
	    					cy = C[lx][ly].getPhantomCompY(m);
	    					// 	ciclo sobre las particulas de la Celda2D vecina
	    					for(nv = 0; nv < C[lxv][lyv].getNP(); ++nv){
	    						ipv = C[lxv][lyv].getPar(nv);
	    						dx = Math.abs(x[ip]-x[ipv]-cx);
	    						dy = Math.abs(y[ip]-y[ipv]-cy);
	    						dr = Math.sqrt((dx*dx)+(dy*dy));
	    						//if((dx < radius) && (dy < radius)){
	    						if(dr < radius){
	    							vx[ip] += Math.cos(a[ipv]);
	    							vy[ip] += Math.sin(a[ipv]);
	    							++nvecinos;
	    						}
	    					}
	    				}
	    				vx[ip] = vx[ip]/nvecinos;
	    				vy[ip] = vy[ip]/nvecinos;
	    				dens += nvecinos;
	    			}
	    		}
	    	}
	    	dens /= N;
	    	for(n = 0; n < N; ++n){
	    		da = 2*Math.PI*Math.random();
	    		dx = vx[n] + eta2*Math.cos(da);
	    		dy = vy[n] + eta2*Math.sin(da);
	    		
	    		a[n] = Math.atan2(dy,dx);
	    		a[n] += (2*Math.PI*eta1*(Math.random()-0.5));
	    		dx = Math.cos(a[n]);
	    		dy = Math.sin(a[n]);
	    		
	    		x[n] = L*Math.random();
	    		y[n] = L*Math.random();
	    	}
	    	return dens;
	    }// Fin de updateVicsekChate
	}
	//------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------
	
	
	public class Celda2D {

	    int nP = 0;
	    int[] par;
	    pXY[] vecino = new pXY[9];
	    pXY[] compens = new pXY[9];

	    public Celda2D(int i, int j, int L){
			int a,b,c,cx,cy;
			c = 0;
			for(int ix = -1; ix <= 1; ++ix){
			    for(int iy = -1; iy <= 1; ++iy){
				a = i + ix;
				b = j + iy;
				cx = 0;
				cy = 0;
				if(a < 0){
				    a += L;
				    cx = -L;
				}
				if(a >= L){
				    a -= L;
				    cx = L;
				}
				if(b < 0){
				    b += L;
				    cy = -L;
				}
				if(b >= L){
				    b -= L;
				    cy = L;
				}
				vecino[c] = new pXY(a,b);
				compens[c] = new pXY(cx,cy);
				++c;
			    }
			}
	    } 
	    public void setNP(int n){
	    	nP = n;
	    	par = new int[n];
	    }
	    public void setPar(int n, int m){
	    	par[n] = m;
	    }
	    public int getNP(){
	    	return nP;
	    }
	    public int getPar(int n){
	    	return par[n];
	    }
	    public int getVecinoX(int n){
	    	return vecino[n].getX();
	    }
	    public int getVecinoY(int n){
	    	return vecino[n].getY();
	    }
	    public int getPhantomCompX(int n){
	    	return compens[n].getX();
	    }
	    public int getPhantomCompY(int n){
	    	return compens[n].getY();
	    }
	}

	/* ***************************************** */


	class pXY
	{
	    int X;
	    int Y;
	    public pXY(int x, int y){
			X = x;
			Y = y;
	    }
	    public int getX(){
	    	return X;
	    }
	    public int getY(){
	    	return Y;
	    }
	}
	
// ----------------------------------
// ----------------------------------

//VARIABLES GLOBALES

	
	JitterMatrix jm;

	int N; // CANTIDAD DE AGENTES EN LA PANTALLA
    double id[];
    double x[];
    double y[];
    double a[];
    double life[];
    double z[];
    double vx [];
    double vy [];
    double vz [];
    double ax [];
    double ay [];
    double az [];
    double mass[];
    double charge[];
	double alpha[];
    double r[];
    double g[];
    double b[];
    double size[];
    int type[];

	int oldN;
    double oldid[];
    double oldx[];
    double oldy[];
    double olda[];
    double oldlife[];
    double oldz[];
    double oldvx [];
    double oldvy [];
    double oldvz [];
    double oldax [];
    double olday [];
    double oldaz [];
    double oldmass[];
    double oldcharge[];
	double oldalpha[];
    double oldr[];
    double oldg[];
    double oldb[];
    double oldsize[];
    int oldtype[];

	float agentColor[] = new float[] {0.f, 0.f, .0f, 0.f};

	float agentSize[];


	double tx[];
    double ty[];

    int L = 10; // Tamaño Circulo inicial???
    double u = 0.00005; // Velocidad Inicial
    double etaVicsek = 0.5;// Ruido Visek inicial
    double etaChate = 0.5; // Ruido Chate inicial
    double radius = 0.4; // Radius inicial
    
    double psi,dens, avgX, avgY, avgA;

 
    Celda2D[][] C = new Celda2D[L][L];
    
    int spX,spY,oX,oY;
 


// Variables de Atributos
	private double velocity = 0.02;
	private double radio= 0.2;
	private double cNoise= 0.2;
	private double vNoise= 0.2;




    // VARIABLES OUTLETS PARA SALIDA DE DATOS
    

    double K ;

// ----------------------------------


	private static final String[] INLET_ASSIST = new String[]{
		"inlet 1 help"
	};
	private static final String[] OUTLET_ASSIST = new String[]{
		"outlet 1 help"
	};
	public Swarms(){

	bail("(mxj vicsekGenerator) you must provide a default”+ ” NUMBER OF AGENTS e.g. [mxj vicsekGenerator 20].");
}
	public Swarms(int _N)
	{

		init(_N);

		declareInlets(new int[]{DataTypes.ALL, DataTypes.INT});
		declareOutlets(new int[]{DataTypes.ALL,DataTypes.ALL});

		declareAttribute("velocity");
		declareAttribute("radio"); 
		declareAttribute("cNoise"); 
		declareAttribute("vNoise"); 

		declareAttribute("agentSize", null, "setAgentSize");
		declareAttribute("agentColor", null, "setAgentColor");	
		declareAttribute("agentType", null, "setAgentType");	

		setInletAssist(INLET_ASSIST);
		setOutletAssist(OUTLET_ASSIST);
	


		JitterMatrix jm = null;

	}
    

    public void init(int _N)
	{

	   	N = _N;

		 id = new double[N];
		 x = new double[N];
		 y = new double[N];
       	 a = new double[N];
		 life = new double[N];
		 z = new double[N];
       	 vx = new double[N];
		 vy = new double[N];
		 vz = new double[N];
       	 ax = new double[N];     
		 ay = new double[N];
		 az = new double[N];
       	 mass = new double[N];
		 charge = new double[N];
		 alpha = new double [N];
    	 r = new double [N];
    	 g = new double [N];
    	 b = new double [N];
     	 size = new double [N];
     	 type = new int [N];

		 oldid = new double[N];
		 oldx = new double[N];
		 oldy = new double[N];
       	 olda = new double[N];
		 oldlife = new double[N];
		 oldz = new double[N];
       	 oldvx = new double[N];
		 oldvy = new double[N];
		 oldvz = new double[N];
       	 oldax = new double[N];     
		 olday = new double[N];
		 oldaz = new double[N];
       	 oldmass = new double[N];
		 oldcharge = new double[N];
		 oldalpha = new double [N];
    	 oldr = new double [N];
    	 oldg = new double [N];
    	 oldb = new double [N];
     	 oldsize = new double [N];
     	 oldtype = new int [N];


		 agentSize  = new float [2];


		 tx = new double[N];
		 ty = new double[N];

		for(int n = 0; n < N; ++n){

		x[n] = Math.random(); // POSICIÓN INICIAL EN X
		y[n] = Math.random(); // POSICIÓN INICIAL EN Y
		a[n] = 2*Math.PI*Math.random(); // A = ?? angulo?
		id[n] = n;
		life[n] =0;
		z[n]= 0;
       	vx[n] =0 ;
		vy[n] =0;
		vz[n] =0;
       	ax[n] =0;     
		ay[n] =0;
		az[n] =0;
       	mass[n] =0;
		charge[n] =0;
		type[n] = 0;

		oldx[n] = Math.random(); // POSICIÓN INICIAL EN X
		oldy[n] = Math.random(); // POSICIÓN INICIAL EN Y
		olda[n] = 2*Math.PI*Math.random(); // A = ?? angulo?
		oldid[n] = n;
		oldlife[n] =0;
		oldz[n]= 0;
       	oldvx[n] =0 ;
		oldvy[n] =0;
		oldvz[n] =0;
       	oldax[n] =0;     
		olday[n] =0;
		oldaz[n] =0;
       	oldmass[n] =0;
		oldcharge[n] =0;
		oldtype[n] = 0;

		// color y tamaño



		alpha[n] =1;
    	r[n] = 1;
    	g[n] = 1;
    	b[n] = 1;
    	size[n] = 1;

		}

		System.arraycopy(x, 0, tx, 0, x.length);
		System.arraycopy(y, 0, ty, 0, y.length);


		dens = ((1.0*N)/(1.0*L*L))*Math.PI*radius*radius;

		// creates cells
		for(int lx = 0; lx < L; ++lx)
		    for(int ly = 0; ly < L; ++ly)
			C[lx][ly] = new Celda2D(lx,ly,L);

	}


    // =====================================

private void generateNewOutputMatrix()
	{
		// ParticleMatrix (necessary planes: ID,life,x,y,z,vx,vy,vz,ax,ay,az,mass,charge)
		// second row have de the last frame info

		if (jm != null)
		jm.freePeer();
		jm = new JitterMatrix(20, "float64", N, 2);
		
		int offset[] = new int[] {0,1};
		
		int dim[] = jm.getDim();

        jm.copyArrayToVectorPlanar(0,0,null,id,dim[0],0);
        jm.copyArrayToVectorPlanar(1,0,null,life,dim[0],0);
        jm.copyArrayToVectorPlanar(2,0,null,tx,dim[0],0);
        jm.copyArrayToVectorPlanar(3,0,null,ty,dim[0],0);
        jm.copyArrayToVectorPlanar(4,0,null,z,dim[0],0);
        jm.copyArrayToVectorPlanar(5,0,null,vx,dim[0],0);
        jm.copyArrayToVectorPlanar(6,0,null,vy,dim[0],0);
        jm.copyArrayToVectorPlanar(7,0,null,vz,dim[0],0);
        jm.copyArrayToVectorPlanar(8,0,null,ax,dim[0],0);
        jm.copyArrayToVectorPlanar(9,0,null,ay,dim[0],0);
        jm.copyArrayToVectorPlanar(10,0,null,az,dim[0],0);
        jm.copyArrayToVectorPlanar(11,0,null,mass,dim[0],0);
        jm.copyArrayToVectorPlanar(12,0,null,charge,dim[0],0);
        jm.copyArrayToVectorPlanar(13,0,null,a,dim[0],0);
        jm.copyArrayToVectorPlanar(14,0,null,alpha,dim[0],0);
        jm.copyArrayToVectorPlanar(15,0,null,r,dim[0],0);
        jm.copyArrayToVectorPlanar(16,0,null,g,dim[0],0);
        jm.copyArrayToVectorPlanar(17,0,null,b,dim[0],0);
        jm.copyArrayToVectorPlanar(18,0,null,size,dim[0],0);
        jm.copyArrayToVectorPlanar(19,0,null,type,dim[0],0);

        jm.copyArrayToVectorPlanar(0,0,offset,oldid,dim[0],0);
        jm.copyArrayToVectorPlanar(1,0,offset,oldlife,dim[0],0);
        jm.copyArrayToVectorPlanar(2,0,offset,oldx,dim[0],0);
        jm.copyArrayToVectorPlanar(3,0,offset,oldy,dim[0],0);
        jm.copyArrayToVectorPlanar(4,0,offset,oldz,dim[0],0);
        jm.copyArrayToVectorPlanar(5,0,offset,oldvx,dim[0],0);
        jm.copyArrayToVectorPlanar(6,0,offset,oldvy,dim[0],0);
        jm.copyArrayToVectorPlanar(7,0,offset,oldvz,dim[0],0);
        jm.copyArrayToVectorPlanar(8,0,offset,oldax,dim[0],0);
        jm.copyArrayToVectorPlanar(9,0,offset,olday,dim[0],0);
        jm.copyArrayToVectorPlanar(10,0,offset,oldaz,dim[0],0);
        jm.copyArrayToVectorPlanar(11,0,offset,oldmass,dim[0],0);
        jm.copyArrayToVectorPlanar(12,0,offset,oldcharge,dim[0],0);
        jm.copyArrayToVectorPlanar(13,0,offset,olda,dim[0],0);
        jm.copyArrayToVectorPlanar(14,0,offset,oldalpha,dim[0],0);
        jm.copyArrayToVectorPlanar(15,0,offset,oldr,dim[0],0);
        jm.copyArrayToVectorPlanar(16,0,offset,oldg,dim[0],0);
        jm.copyArrayToVectorPlanar(17,0,offset,oldb,dim[0],0);
        jm.copyArrayToVectorPlanar(18,0,offset,oldsize,dim[0],0);
        jm.copyArrayToVectorPlanar(19,0,offset,oldtype,dim[0],0);

		outlet(0, "jit_matrix", jm.getName());
		
		System.arraycopy(id, 0, oldid, 0, id.length);
		System.arraycopy(life, 0, oldlife, 0, life.length);
		System.arraycopy(tx, 0, oldx, 0, x.length);
		System.arraycopy(ty, 0, oldy, 0, y.length);
		System.arraycopy(z, 0, oldz, 0, z.length);
		System.arraycopy(vx, 0, oldvx, 0, vx.length);
		System.arraycopy(vy, 0, oldvy, 0, vy.length);
		System.arraycopy(vz, 0, oldvz, 0, vz.length);
		System.arraycopy(ax, 0, oldax, 0, ax.length);
		System.arraycopy(ay, 0, olday, 0, ay.length);
		System.arraycopy(az, 0, oldaz, 0, az.length);
		System.arraycopy(mass, 0, oldmass, 0, mass.length);
		System.arraycopy(charge, 0, oldcharge, 0, charge.length);

		System.arraycopy(a, 0, olda, 0, a.length);

		System.arraycopy(alpha, 0, oldalpha, 0, alpha.length);
		System.arraycopy(r, 0, oldr, 0, r.length);
		System.arraycopy(g, 0, oldg, 0, g.length);
		System.arraycopy(b, 0, oldb, 0, b.length);

		System.arraycopy(type, 0, oldtype, 0, type.length);


	} 

private void calculateAverage()
	{
		double sumX=0;
		double sumY=0;
		double sumA=0;
		for(int n = 0; n < N; ++n)
		{

		// Escalado de (0 - 10) - (-1 _ +1) 
		tx[n] = (x[n]*0.2)-1;
		ty[n] = (y[n]*0.2)-1;

		sumX = sumX + tx[n];
		sumY = sumY + ty[n];
		sumA = sumA + a[n];

		}
		
		avgX = sumX / x.length;
		avgY = sumY / y.length;
		avgA = sumA / a.length;
	}

	private void setAgentColor(float[] _agentColor)
	{




	if (_agentColor.length == 5)
	{
	int agentChange = (int)_agentColor[0]; // agentID
			if (agentChange >= 0 && agentChange < N )// limites de id
				{
				r[agentChange]= _agentColor[1];
				g[agentChange]= _agentColor[2];
				b[agentChange]= _agentColor[3];
				alpha[agentChange]= _agentColor[4];

			}
		
	}
	if (_agentColor.length == 4)
			{
		for(int i=0; i<r.length; i++)
				{

				r[i]= _agentColor[0];
				g[i]= _agentColor[1];
				b[i]= _agentColor[2];
				alpha[i]= _agentColor[3];

				}
			}	
	}

	public void randomColor()
	{



	Random randC = new Random();

		for(int i=0; i<r.length; i++)
				{

				r[i]= (float)randC.nextGaussian();
				g[i]= (float)randC.nextGaussian();
				b[i]= (float)randC.nextGaussian();
				alpha[i]= (float)randC.nextGaussian();
				


				}

	}


	private void setAgentSize( float[] _agentSize)
	{


	if (_agentSize.length == 2)
		{

		if (_agentSize[1] > 0 )// limites de valor
			{

			if (_agentSize[0] >= 0 && _agentSize[0] < N ) // limites de id
				{
				size [(int)_agentSize[0]] =  _agentSize[1];
			
				System.arraycopy(size, 0, oldsize, 0, size.length);

				}
			}
		}


	if (_agentSize.length == 1)
		{

		if (_agentSize[0] > 0 )// limites de valor	
			{


		for(int i=0; i<size.length; i++)
   			 	{
		size[i]= _agentSize[0];
//		post ("estoy en cambiando TODOS" +"_ _ "+  size [0]);
		System.arraycopy(size, 0, oldsize, 0, size.length);

				}

		
			}
		}

	}

	public void randomSize()
	{



	Random randS = new Random();

		for(int i=0; i<r.length; i++)
				{

				size[i]= 1.f + Math.abs((float)randS.nextGaussian())*30;


				}

	}
	public void setAgentType(float[] _agentType)
	{

	if (_agentType.length == 2)
		{

		if (_agentType[1] >= 0 )// limites de valor
			{

			if (_agentType[0] >= 0 && _agentType[0] < N ) // limites de id
				{
				type [(int)_agentType[0]] =  (int)_agentType[1];
			
				}
			}
		}


	if (_agentType.length == 1)
		{

		if (_agentType[0] >= 0 )// limites de valor	
			{


		for(int i=0; i<type.length; i++)
   			 	{
		type[i]= (int)_agentType[0];

				}

		
			}
		}

	
	}
	public void randomType()
	{



	Random randT = new Random();

		for(int i=0; i<r.length; i++)
				{

				type[i]= (int)(Math.abs(randT.nextGaussian())*3);


				}

	}


	public void getNumberOfParticles()
	{
	outlet (1, new Atom[]{Atom.newAtom("numberOfParticles/"), Atom.newAtom(N)});
	}


	public void bang()
	{

	u = velocity;
	radius = radio;
	etaVicsek= vNoise;
	etaChate = cNoise;
	
	generateNewOutputMatrix();

	Methods2D.updateCeldas(C,x,y);
	dens = Methods2D.updateVicsekChate(C,x,y,a,radius,etaVicsek,etaChate,u);
	psi = Methods2D.updatePsi(a);
	calculateAverage();

	outlet (1, new Atom[]{Atom.newAtom(psi), Atom.newAtom(dens),Atom.newAtom(avgX),Atom.newAtom(avgY),Atom.newAtom(avgA)});
	
	oldN = N; 

	}

	

    
	public void inlet(int i)
	{
		if (getInlet() == 1) {
		init(i);
		}
	}
    
	public void inlet(float f)
	{
	
	}
    
    
	public void list(Atom[] list)
	{
	}
    
}

    // =====================================





