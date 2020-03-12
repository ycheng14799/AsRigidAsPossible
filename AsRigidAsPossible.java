import java.awt.*;
import java.awt.event.*;
import java.applet.Applet;
import java.io.*;

/* Class AsRigidAsPossible */ 
public class AsRigidAsPossible extends Applet {
	static Frame myFrame = null; 
	Button bTriangulate; 
	Checkbox cAnimate; 
	Button bClear; 
	Button bQuit; 
	Panel mainPanel; 
	private MyGraphics myG; 
	public int borderSize; 
	public Color paintColor, bkColor, constraintColor; 

	public AsRigidAsPossible() { // Constructor
		borderSize = 10000;
        paintColor = Color.red;
        bkColor    = Color.lightGray;
		constraintColor = Color.blue;
		bTriangulate = new Button("Triangulate");
		bClear = new Button("Clear");
		cAnimate = new Checkbox("Animate");
	}

	public void clearMe() {
		Graphics g = getGraphics(); 
		Dimension dimension = getSize();
		g.setColor(getBackground());
		g.fillRect(0, 0, dimension.width, dimension.height);
	}

	public void init() { // Initialize the applet 
		setBackground(Color.gray);
		setLayout(new BorderLayout());
		
		// Main panel 
		mainPanel = new Panel(); 
			mainPanel.setBackground(Color.lightGray);
		
		// UI
		Panel panel2 = new Panel(); 
			panel2.setBackground(Color.lightGray);
		bTriangulate = new Button("Triangulate");
			panel2.add(bTriangulate);
		cAnimate = new Checkbox("Animate");
			panel2.add(cAnimate);
		bClear = new Button("Clear"); 
			panel2.add(bClear);
		if (myFrame != null) {
			bQuit = new Button("Quit");
			panel2.add(bQuit);
		}

		add("North", panel2);
		add("South", mainPanel);

		myG = new MyGraphics(paintColor, Color.yellow, Color.gray);
		
		// Triangulation function
		bTriangulate.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				clearMe();
				Graphics g = getGraphics();
				myG.triangulate(g);
				g.setPaintMode();
			}
		});

		// Prepare for animation 
		cAnimate.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {               
				if(e.getStateChange() == 1) {
					myG.preComputeForManipulation(); 
				}
			} 
		});

		// Clear function 
		bClear.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				// Clear components 
				myG.clear(mainPanel, getBackground()); 
				// Clear screen 
				clearMe(); 
				Graphics g = getGraphics(); 
				g.setPaintMode();
			}
		});

		// Quit function 
		if (myFrame != null) {
		  bQuit.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				   System.exit(0);
			}
		  });
		}

		// Handling mouse input 
		this.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				int x = e.getX();
				int y = e.getY();
				if(!myG.thePoly.triangulated) {
					Graphics g = getGraphics(); 
					g.setColor(paintColor);
					myG.addVertex(g, x, y);
					g.setPaintMode(); 
				} else if(!cAnimate.getState()) {
					Graphics g = getGraphics(); 
					g.setColor(constraintColor);
					myG.addConstraint(g, x, y);
					g.setPaintMode();
				}
			}
		});


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