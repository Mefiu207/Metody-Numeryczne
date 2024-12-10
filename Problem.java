import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;


public class Problem extends Thread{
    
    public static Integer IT_MAX = 500;
    
    private Double D;

    private Double dt;
    private Double delta;

    private Integer nx;
    private Integer ny;

    private Double[] x;
    private Double[] y;

    private Integer i1;
    private Integer i2;
    private Integer j1;

    private Double sigma;

    private Double xa;
    private Double ya;

    private Double psi[][];

    private Double vx[][];
    private Double vy[][];

    private Double u0[][];
    private Double u1[][];

    private ArrayList<Double> c;
    private ArrayList<Double> xsr;



    public Problem(Double Dval){
        D = Dval;
        
        delta = 0.01;

        nx = 400;
        ny = 90;

        i1 = 200;
        i2 = 210;
        j1 = 50;

        sigma = 10 * delta;

        xa = 0.45;
        ya = 0.45;

        // Initialize arrays
        x = new Double[nx + 1];
        y = new Double[ny + 1];

        psi = new Double[ny + 1][nx + 1];
        vx = new Double[ny + 1][nx + 1];
        vy = new Double[ny + 1][nx + 1];

        u0 = new Double[ny + 1][nx + 1];
        u1 = new Double[ny + 1][nx + 1];

        c = new ArrayList<Double>();
        xsr = new ArrayList<Double>();

        populateXandY();
    }

    // Populates x and y arrays with correct values
    private void populateXandY(){
        for(int i = 0; i <= nx; i++) x[i] = i * delta;
        for(int j = 0; j <= ny; j++) y[j] = j * delta;
    }

    // Reades data from a given file to psi array
    public void readPsiFromFile(String filename){
        
        try{    
            File myFile = new File(filename);
            Scanner myScanner = new Scanner(myFile);

            while(myScanner.hasNextLine()){

                String line = myScanner.nextLine().trim();

                String[] parts = line.split("\\s+");

                Integer xPoint = Integer.parseInt(parts[0]);
                Integer yPoint = Integer.parseInt(parts[1]);
                
                Double value = Double.parseDouble(parts[2]);

                psi[yPoint][xPoint] = value;
            }
            myScanner.close();
        
        } catch(FileNotFoundException e) { e.printStackTrace(); }

    }

    // Calculates values for Vx and Vy arrays
    public void calcV(){
        
        // Inner field
        for(int j = 1; j < ny; j++){
            for(int i = 1; i < nx; i++){
                
                vx[j][i] = (psi[j + 1][i] - psi[j - 1][i]) / (2 * delta);
                vy[j][i] = -(psi[j][i + 1] - psi[j][i - 1]) / (2 * delta);

            }
        }

        // Edges inside the box
        for(int j = 0; j <= j1; j++){
            for(int i = i1; i <= i2; i++){

                vx[j][i] = 0.0;
                vy[j][i] = 0.0;

            }
        }

        // Upper and lower edge
        for(int i = 1; i < nx; i++){

            vx[0][i] = 0.0;
            vy[0][i] = 0.0;

            vx[ny][i] = 0.0;
            vy[ny][i] = 0.0;

        }

        // Left and right edge
        for(int j = 0; j <= ny; j++){

            vx[j][0] = vx[j][1];
            vx[j][nx] = vx[j][nx - 1];

            vy[j][0] = vy[j][1];
            vy[j][nx] = vy[j][nx - 1];

        }

        setDeltaT();

    }    

   // Finds Vmax and sets value of dt (delta T) (calulated in calV() method)
    private void setDeltaT(){
        Double currentHighest = 0.0;

        Double presentV = 0.0;

        for(int j = 1; j < ny; j++){
            for(int i = 1; i < nx; i++){

                presentV = Math.sqrt(Math.pow(vx[j][i], 2) + Math.pow(vy[j][i], 2)); 
                if(presentV > currentHighest) currentHighest = presentV;
                
            }
        }

        this.dt = delta/(4.0*currentHighest);
    }

    // Saves values of vx and vy to given file
    public void saveVToFile(String filename){
        
        File myFile = new File(filename);
        myFile.delete();

        try {myFile.createNewFile();}
        catch(IOException e) {e.printStackTrace();}

        try{
            FileWriter writeFile = new FileWriter(filename);

            for(int j = 0; j <= ny; j++){
                for(int i = 0; i <= nx; i++){
                    writeFile.write(i + " " + j + " " + vx[j][i] + " " + vy[j][i] + "\n");
                }
                writeFile.write("\n");
            }

            writeFile.close();
        }catch(IOException e) { e.printStackTrace(); }   
    } 

    // Calculates the first u0 value 
    private void setU0(){

        Double firstPart = 1.0/(2.0*Math.PI*Math.pow(sigma, 2));

        for(int j = 0; j <= ny; j++){
            for(int i = 0; i <= nx; i++){

                u0[j][i] = firstPart * Math.exp(-(Math.pow(x[i] - xa, 2) + Math.pow(y[j] - ya, 2)) / (2.0 * Math.pow(sigma, 2))); 

            }
        }

    }

    // Cheks if given indexes are part of the wall in the midle of the net 
    private boolean checkIfEdge(int i, int j){
        if(i >= i1 && i <= i2 && j >= 0 && j <= j1) return true;
        else return false;
    }

    // Main iterating function
    private void picardIter(){
        setU0();

        Integer kTab[] = {1, 2, 3, 4, 5};

        Double Tk[] = new Double[5];

        


        for(int IT = 1; IT <= IT_MAX; IT++){
            u1 = u0.clone();

            for(int k = 1; k <= 20; k++){
                for(int i = 0; i <= nx; i++){
                    for(int j = 1; j <= ny - 1; j++){
                        
                        if(checkIfEdge(i, j)) continue;
                        else calcNextU(i, j);

                    }
                }
            }
            u0 = u1.clone();
            
            calcC();
            calcXsr();
        



        }

    }

    // Checks if periodic edges are in place and calculates a value for u1 of given indexes
    private void calcNextU(int i, int j){
        
        Double frontValue = 1.0 / (1.0 + (2.0*D * dt / Math.pow(delta, 2)));
        
        Double dt_vx = 0.0;
        Double dt_vy = 0.0;
        Double dt_D = 0.0;

        if(i == 0){
            dt_vx = dt/2.0* vx[j][i] * (((u0[j][i + 1] - u0[j][nx]) / (2.0*delta)) + ((u1[j][i + 1] - u1[j][nx]) / (2.0*delta)));
        
            dt_vy = dt/2.0* vy[j][i] * (((u0[j + 1][i] - u0[j - 1][i]) / (2.0*delta)) + ((u1[j + 1][i] - u1[j - 1][i]) / (2.0*delta)));
        
            dt_D = dt/2.0* D * ( (u0[j][i + 1] + u0[j][nx] + u0[j + 1][i] + u0[j - 1][i] - 4*u0[j][i]) / Math.pow(delta, 2) + (u1[j][i + 1] + u1[j][nx] + u1[j + 1][i] + u1[j - 1][i]) / Math.pow(delta, 2));
        
        }else if(i == nx){
            dt_vx = dt/2.0* vx[j][i] * (((u0[j][0] - u0[j][i - 1]) / (2.0*delta)) + ((u1[j][0] - u1[j][i - 1]) / (2.0*delta)));
        
            dt_vy = dt/2.0* vy[j][i] * (((u0[j + 1][i] - u0[j - 1][i]) / (2.0*delta)) + ((u1[j + 1][i] - u1[j - 1][i]) / (2.0*delta)));
        
            dt_D = dt/2.0* D * ( (u0[j][0] + u0[j][i - 1] + u0[j + 1][i] + u0[j - 1][i] - 4*u0[j][i]) / Math.pow(delta, 2) + (u1[j][0] + u1[j][i - 1] + u1[j + 1][i] + u1[j - 1][i]) / Math.pow(delta, 2));

        }else{
            dt_vx = dt/2.0* vx[j][i] * (((u0[j][i + 1] - u0[j][i - 1]) / (2.0*delta)) + ((u1[j][i + 1] - u1[j][i - 1]) / (2.0*delta)));
        
            dt_vy = dt/2.0* vy[j][i] * (((u0[j + 1][i] - u0[j - 1][i]) / (2.0*delta)) + ((u1[j + 1][i] - u1[j - 1][i]) / (2.0*delta)));
        
            dt_D = dt/2.0* D * ( (u0[j][i + 1] + u0[j][i - 1] + u0[j + 1][i] + u0[j - 1][i] - 4*u0[j][i]) / Math.pow(delta, 2) + (u1[j][i + 1] + u1[j][i - 1] + u1[j + 1][i] + u1[j - 1][i]) / Math.pow(delta, 2));

        }

        u1[j][i] = frontValue * ( u0[j][i] - dt_vx - dt_vy + dt_D);
    }

    // Calculates current value of c
    private void calcC(){
        
        Double value = 0.0;
        
        for(int j = 0; j <= ny; j++){
            for(int i = 0; i <= nx; i++){
                value += u0[j][i] * Math.pow(delta, 2);
            }
        }

        c.add(value);
    }

    // Calculates current value of xsr
    private void calcXsr(){

        Double value = 0.0;
        
        for(int j = 0; j <= ny; j++){
            for(int i = 0; i <= nx; i++){
                value += u0[j][i] * Math.pow(delta, 2) * x[i];
            }
        }

        xsr.add(value);

    }

    // Saves C data to given file
    public void saveCandXSRtoFile(String filename){
           
        File myFile = new File(filename);
        myFile.delete();

        try {myFile.createNewFile();}
        catch(IOException e) {e.printStackTrace();}

        try{
            FileWriter writeFile = new FileWriter(filename);

            Double Ndt = 0.0;

            for(int j = 0; j < c.size(); j++){
                
                Ndt = j*dt;
                writeFile.write( Ndt + " " + c.get(j) + " " + xsr.get(j) + "\n");
            }

            writeFile.close();
        }catch(IOException e) { e.printStackTrace(); }   

    }

    // Function for using threads that starts picardIter function
    public void run(){
        picardIter();
    }

    public static void main(String[] args){
        
        long begin = System.nanoTime();
        long now;

        Problem a1 = new Problem(0.0);
        
        a1.readPsiFromFile("psi.dat");
        a1.calcV();
        a1.saveVToFile("v.dat");

        System.out.println("First part is starting...");


        Problem a2 = new Problem(0.1);
        a2.readPsiFromFile("psi.dat");
        a2.calcV();

        System.out.println("Second part is starting...");

        a1.start();
        System.out.println("Starting first thred");

        a2.start();
        System.out.println("Starting second thred");

        while(a1.isAlive() || a2.isAlive()){

        }

        a1.saveCandXSRtoFile("ct_d0.dat");
        a2.saveCandXSRtoFile("ct_d01.dat");

        now = System.nanoTime();
        
        long elapsedTimeNano = now - begin;
        long elapsedTimeSeconds = elapsedTimeNano / 1_000_000_000;
        long minutes = elapsedTimeSeconds / 60;
        long seconds = elapsedTimeSeconds % 60;
        
        System.out.println("Time the program took with " + IT_MAX + " iterations:");
        System.out.println(minutes + " m " + seconds + " s");

    }
}