package Swarms;

import com.cycling74.max.*;
import com.cycling74.jitter.*;

public class Swarms2data extends MaxObject
{
	JitterMatrix output = new JitterMatrix();
	
	public int[] dim;
    double id[]; 
	double x[];
    double y[];
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
    double a[];

	public int numberOfParticles = 0;


	private static final String[] OUTLET_ASSIST = new String[]{
		"outlet 1 help"
	};
	public Swarms2data(float gain)
	{
		declareIO(1,1);
		setInletAssist(0, "input (matrix)");
		setOutletAssist(OUTLET_ASSIST);
	}
    
	public void jit_matrix(String mname)
	{
		output.frommatrix(mname);

		dim = output.getDim();
		numberOfParticles = dim[0];




        id = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(0,0,null,x,dim[0],0);
        life = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(1,0,null,x,dim[0],0);
        x = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(2,0,null,x,dim[0],0);
        y = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(3,0,null,y,dim[0],0);
        z = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(4,0,null,z,dim[0],0);
        vx = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(5,0,null,vx,dim[0],0);
        vy = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(6,0,null,vy,dim[0],0);
        vz = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(7,0,null,vz,dim[0],0);
        ax = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(8,0,null,ax,dim[0],0);
        ay = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(9,0,null,ay,dim[0],0);
        az = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(10,0,null,az,dim[0],0);
        mass = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(11,0,null,mass,dim[0],0);
        charge = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(12,0,null,charge,dim[0],0);
		a = new double[numberOfParticles];
        output.copyVectorToArrayPlanar(13,0,null,a,dim[0],0);

		for (int i=0;i<dim[0];i++){
		outletHigh(0,new Atom[]{Atom.newAtom(i),Atom.newAtom( x[i]),Atom.newAtom( y[i]) ,Atom.newAtom( a[i])});
			}
		}

}







