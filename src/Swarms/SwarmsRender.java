package Swarms;

import com.cycling74.max.*;
import com.cycling74.jitter.*;


public class SwarmsRender extends MaxObject
{

	JitterMatrix input = new JitterMatrix();
	JitterMatrix output = new JitterMatrix();

	JitterObject window;
	JitterObject render;
	JitterObject sketch;


	int shape = 0, drawMode = 1, L, oldL;
	
	

	float backgroundColor[] = new float[] {1.f, 1.f, 1.f, .1f};






	public SwarmsRender()
		{


		declareAttribute("shape");
		declareAttribute("drawMode");
		declareAttribute("backgroundColor");


		declareInlets(new int[]{DataTypes.ALL});
		declareOutlets(new int[]{DataTypes.ALL});
   
		init();



		}

	public void jit_matrix(String mname)
	{
		input.frommatrix(mname);	
		int dim[] = input.getDim();
		L = dim[0];

		double x[] = new double[dim[0]];
		double y[] = new double[dim[0]];
		double a[] = new double[dim[0]];
		double alpha[] = new double[dim[0]];
		double r[] = new double[dim[0]];
		double g[] = new double[dim[0]];
		double b[] = new double[dim[0]];
		float size[] = new float[dim[0]];


		//oldargb = new float [dim[0]][4];				

		input.copyVectorToArrayPlanar(2,0,null,x,dim[0],0);
        input.copyVectorToArrayPlanar(3,0,null,y,dim[0],0);
		input.copyVectorToArrayPlanar(13,0,null,a,dim[0],0);
		input.copyVectorToArrayPlanar(14,0,null,alpha,dim[0],0);
        input.copyVectorToArrayPlanar(15,0,null,r,dim[0],0);
		input.copyVectorToArrayPlanar(16,0,null,g,dim[0],0);
        input.copyVectorToArrayPlanar(17,0,null,b,dim[0],0);
        input.copyVectorToArrayPlanar(18,0,null,size,dim[0],0);				

	
		render.setAttr("erase_color",backgroundColor);

		sketch.send("reset");





		for(int n = 0; n < dim[0]; ++n)
		{
		
	
		sketch.send("glpointsize",new float[]{size[n]});	

		sketch.send("glcolor",new float[]{(float)r[n],(float)g[n],(float)b[n]});			

		sketch.send("point", new float[]{(float)x[n],(float)y[n],(float)0.});

		}

		oldL = L;
	}

	public void bang()
    {

		render.send("erase");
		render.send("drawclients");
		render.send("swap");
    }




	public void init()
	{

		// create our window
		window = new JitterObject("jit.window",new Atom[]{Atom.newAtom("swarmsWindow")});
		window.setAttr("size", new int[] {640, 480} ); 	
		window.setAttr("depthbuffer",1); 
		window.setAttr("fsaa",1); //full scene anti-alliasing

		// create our render object for our window
		render = new JitterObject("jit.gl.render",new Atom[]{Atom.newAtom("swarmsWindow")});
		render.setAttr("camera",new float[] {0.f,0.f,4.f});

		// create out sketch object
		sketch = new JitterObject("jit.gl.sketch",new Atom[]{Atom.newAtom("swarmsWindow")});
		sketch.setAttr("lighting_enable",1);
		sketch.setAttr("smooth_shading",0);
		sketch.setAttr("displaylist",1);
	
		

	}
	public void notifyDeleted()
	{
		sketch.freePeer();
		render.freePeer();
		window.freePeer();
	}
    
}

