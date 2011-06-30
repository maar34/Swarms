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


	int shape = 0, agentSize = 1, drawMode = 1;

	float backgroundColor[] = new float[] {1.f, 1.f, 1.f, .1f};
	float agentColor[] = new float[] {0.f, 0.f, .0f, 0.f};


	public SwarmsRender()
		{


		declareAttribute("shape");
		declareAttribute("agentSize");
		declareAttribute("drawMode");
		declareAttribute("backgroundColor");
		declareAttribute("agentColor");

		declareInlets(new int[]{DataTypes.ALL});
		declareOutlets(new int[]{DataTypes.ALL});
   
		init();



		}

	public void jit_matrix(String mname)
	{
		input.frommatrix(mname);	
		int dim[] = input.getDim();

		double x[] = new double[dim[0]];
		double y[] = new double[dim[0]];
		double a[] = new double[dim[0]];

		input.copyVectorToArrayPlanar(2,0,null,x,dim[0],0);
        input.copyVectorToArrayPlanar(3,0,null,y,dim[0],0);
		
		render.setAttr("erase_color",backgroundColor);

		sketch.send("reset");

	//	sketch.send("glPointSize",(float)agentSize);			
		sketch.send("glcolor",agentColor);	
		for(int n = 0; n < dim[0]; ++n)
		{
		
		sketch.send("point", new float[]{(float)x[n],(float)y[n],(float)0.});
		}


		//render.send("erase");
		//render.send("drawclients");
		//render.send("swap");

	}

	public void bang()
    {
		sketch.send("glcolor", agentColor);
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