import com.cycling74.max.*;
import com.cycling74.jitter.*;
import java.awt.*;


public class SwarmsRender extends MaxObject
{

	JitterMatrix input = new JitterMatrix();
	JitterMatrix output = new JitterMatrix();

	JitterObject window;
	JitterObject render;
	JitterObject sketch;


	int shape = 0, drawMode = 1, L, oldL;
	
	

	float backgroundColor[] = new float[] {1.f, 1.f, 1.f, .1f};

	float agentColor[] = new float[] {0.f, 0.f, .0f, 0.f};

	float argb [][]; //oldargb[][]

	float agentSize[];

	float agentSizeArray[][];



	public SwarmsRender()
		{


		declareAttribute("shape");
		declareAttribute("agentSize");
		declareAttribute("drawMode");
		declareAttribute("backgroundColor");
		declareAttribute("agentColor");//, null, "setAgentColor"

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

		//oldargb = new float [dim[0]][4];				

		input.copyVectorToArrayPlanar(2,0,null,x,dim[0],0);
        input.copyVectorToArrayPlanar(3,0,null,y,dim[0],0);
		

		if (L!=oldL)
			{
			argb = new float [L][4];
			agentSizeArray = new float [L][1];
			agentSize  = new float [2];

			for(int n = 0; n < L; ++n) 
				{	
				agentSizeArray [n][0] = 1.f;
				}
			}

		render.setAttr("erase_color",backgroundColor);

		sketch.send("reset");



		setColor (agentColor);
	
		setSize (agentSize);

//		sketch.send("glpointsize",(float)agentSize);	

		for(int n = 0; n < dim[0]; ++n)
		{
		
		sketch.send("glpointsize",agentSizeArray[n][0]);	

		sketch.send("glcolor",argb[n]);	

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


	private void setColor(float[] _agentColor)
	{
	int agentChange = (int)_agentColor[0]; // agentID



	if (_agentColor.length == 5)
	{

			if (agentChange >= 0 && agentChange < L )// limites de id
				{
				for(int k=1; k<5; k++)
				
    			{
				argb[agentChange][k-1]= _agentColor[k];
				}
			}
		
	}
	if (_agentColor.length == 4)
			{
	for(int i=0; i<argb.length; i++)
				{
  		for(int j=0; j<argb[i].length; j++)
   		 			{
		
		argb[i][j]= _agentColor[j];


					}
				}
			}
	
	}

	private void setSize( float[] _agentSize)
	{



	if (_agentSize.length == 2)
		{
		if (_agentSize[1] > 0 )// limites de valor
			{
			if (_agentSize[0] >= 0 && _agentSize[0] < L ) // limites de id
				{
				agentSizeArray[(int)_agentSize[0]][0] =  _agentSize[1];
				}
			}
		}

	if (_agentSize.length == 1)
		{

	if (_agentSize[0] > 0 )// limites de valor	
		{
		for(int i=0; i<agentSizeArray.length; i++)
   		 	{
		
		agentSizeArray[i][0]= _agentSize[0];
			}


		
		}
	}

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

