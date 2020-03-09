import java.awt.*;
import java.awt.event.*;
import java.applet.Applet;
import java.io.*;

/* Class AsRigidAsPossible */ 
public class AsRigidAsPossible extends Applet {
	static Frame myFrame = null; 
	Button bQuit; 
	Panel mainPanel; 

	public AsRigidAsPossible() { // Constructor
	
	}

	public void init() { // Initialize the applet 
		
	}
	
	public static void main(String[] args) { // Main
		myFrame = new Frame("As-Rigid-As-Possible Shape Manipulation Application"); // Create frame with title
		AsRigidAsPossible myApp = new AsRigidAsPossible(); // Initialize application
		myApp.init(); 

		// Add applet to the frame
		myFrame.add(myApp, BorderLayout.CENTER); // Add applet to frame 
		myFrame.pack(); // Set window size 
		myFrame.setSize(600, 500);
		myFrame.setVisible(true); // Make frame visible 

	} 
}

/*
public class Bezier extends Applet {
    
    public static final long serialVersionUID = 24362462L;
    static Frame myFrame=null;
    Button bClear;
    Button bRead;
    Button bElevate;
	Button bElevateMult;
	TextField bElevateEntry; 
    Button bQuit;
    Panel  mainPanel;
    private MyGraphics myG;
    public Color paintColor, bkColor;
    public int borderSize;
	private boolean pointSelected; 
	private int selectedIndex; 
    
    
    public Bezier()
    {
        borderSize = 10000;
        paintColor = Color.red;
        bkColor    = Color.lightGray;
        bClear     = new Button("Clear all");
        bRead      = new Button("Read sample points");
        bElevate    = new Button("Elevate once");
		bElevateMult = new Button("Elevate multiple");
		bElevateEntry = new TextField("Num. Elevate");
    }

    public void clearMe() {
	Graphics g = getGraphics();
        Dimension dimension = getSize();
        Color col = g.getColor();
        g.setColor(getBackground());
        g.fillRect(0, 0, dimension.width, dimension.height);
        g.setColor(col);
    }

    public void init() {
	//	System.out.println("Init method");

	setBackground(Color.gray);
        setLayout(new BorderLayout());
        mainPanel  = new Panel();
	mainPanel.setBackground(Color.lightGray);
        Panel panel2 = new Panel();
	panel2.setBackground(Color.lightGray);
        bRead = new Button("Read sample points");
	panel2.add(bRead);
        bElevate = new Button("Elevate once");
	panel2.add(bElevate);
		bElevateEntry = new TextField("Num. Elevate"); 
	panel2.add(bElevateEntry);
		bElevateMult = new Button("Elevate multiple");
	panel2.add(bElevateMult);
        bClear = new Button("Clear all");
	panel2.add(bClear);
    if (myFrame != null)
    {
           bQuit = new Button("Quit");
           panel2.add(bQuit);
    }

	add("North",  panel2);
	add("South",  mainPanel);
        myG = new MyGraphics(paintColor, Color.yellow);
    

	bRead.addMouseListener(new MouseAdapter() {
	    public void mouseClicked(MouseEvent e) {
		int x, y;

		File f = new File(".","sample.data");
		if (!f.exists()) throw new Error("Sample Point File Not Found"); 
		else {
		    // get rid of existing polygon
		    myG.clear(mainPanel, getBackground());
		    // clear screen
		    Graphics g = getGraphics();
		    g.setColor(paintColor);
		    clearMe();
		    // read in points from file
		    try{
			Reader r = new BufferedReader(
				   new InputStreamReader(
				   new FileInputStream(f)));
			StreamTokenizer st = new StreamTokenizer(r);
			
			int i=0;
			while (st.nextToken()!=st.TT_EOF) {
			    x = (int)(st.nval);
			    st.nextToken();
			    y = (int)(st.nval);
			    System.out.println("Read point (" + x + " , "+y+")");
			    myG.addPolyPoint(g, x, y);
			}		
		    }
		    catch(Exception exc){
			System.out.println("Cannot read from file"+f);
		    }
		    g.setPaintMode();
		}
	    }
	});
	
	bElevate.addMouseListener(new MouseAdapter() {
	    public void mouseClicked(MouseEvent e) {
		Graphics g = getGraphics();
		g.setColor(paintColor);
		clearMe();
		myG.elevateOnce(g);
		g.setPaintMode();
	  }
	});//elevate once

	bElevateMult.addMouseListener(new MouseAdapter() {
		public void mouseClicked(MouseEvent e) {
			Graphics g = getGraphics();
			g.setColor(paintColor);
			clearMe();
			String numElevString = bElevateEntry.getText();
			int numElev = Integer.parseInt(numElevString);
			myG.elevateMultiple(g, numElev);
			System.out.println("Elevating " + numElev + " times");
			g.setPaintMode();
		}
	}); // elevate multiple times

	bClear.addMouseListener(new MouseAdapter() {
	    public void mouseClicked(MouseEvent e) {
		// get rid of components
		myG.clear(mainPanel, getBackground());
		// clear screen
		Graphics g = getGraphics();
		g.setColor(paintColor);
		//System.out.println("Clear all");
		clearMe();
		g.setPaintMode();
	    }
	});
    if (myFrame != null)
    {//Quit button
      bQuit.addMouseListener(new MouseAdapter() {
        public void mouseClicked(MouseEvent e) {
               System.exit(0);
        }
      });
    }


	this.addMouseListener(new MouseAdapter() {
		public void mousePressed(MouseEvent e) {
			Graphics g = getGraphics();
			g.setColor(paintColor);
			int x = e.getX();
			int y = e.getY();
			int[] currPointsX = myG.thePoly.xpoints;
			int[] currPointsY = myG.thePoly.ypoints;
			int numCurrPoints = myG.thePoly.npoints; 
			pointSelected = false;
			double pointDist = 0; 
			for(int i=0; i<numCurrPoints;i++) {
				pointDist = (x-currPointsX[i])*(x-currPointsX[i])
					+ (y-currPointsY[i])*(y-currPointsY[i]);
				System.out.println(pointDist);
				if(pointDist < 30.0f) {
					pointSelected = true; 
					selectedIndex = i;
					break; 
				}
			}
			if (pointSelected  == false) {
				myG.addPolyPoint(g, x, y);
				g.setPaintMode();
			} else {
				System.out.println("Point selected");
			}
		}
		public void mouseReleased(MouseEvent e) {
			if (pointSelected) { 
				System.out.println("Dropping point");
				int x = e.getX();
				int y = e.getY(); 
				clearMe();
				Graphics g = getGraphics();
				myG.updatePoint(g, selectedIndex, x, y);
				g.setPaintMode();
			}
		}
	});
    }
	

  public void update(Graphics g) { 
    paint(g); 
  }
  
  public void paint(Graphics g) {
    paintComponents(g); 
    //    System.out.println("MAIN PAINT");
    myG.redrawThePolygon(g);
  }

  public void stop() {
    // Stop the bouncer thread if necessary.
      System.out.println("stop");
  }  

  public static void main(String[] args) {
    Bezier myBezierApplet = new Bezier(); 
    myFrame = new Frame("Bezier degree application"); // create frame with title
    // Call applet's init method (since Java App does not
    // call it as a browser automatically does)
    myBezierApplet.init(); 

    // add applet to the frame
    myFrame.add(myBezierApplet, BorderLayout.CENTER);
    myFrame.pack(); // set window to appropriate size (for its elements)
    myFrame.setSize(600, 500);
    myFrame.setVisible(true); // usual step to make frame visible

  } // end main

}
*/