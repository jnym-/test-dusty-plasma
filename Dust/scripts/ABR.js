/*jshint esversion: 6 */
var dom1 = {
    gSlider: $("input#gammaSec1"),//gamma slider
    JSlider: $("input#jnormSec1"),//j slider
};
var plt1 = {//layout of graph
        layoutJ: {
            autosize: true,
            showlegend: false,
            xaxis: {
                range: [0, 0.5],
                title: "$\\Phi_b$",
            },
            yaxis: {
                type: 'log',
                title: "$J/\\gamma$",
                range: [-10, 10],
            },
            margin: {
                l: 50, r: 10, b:50, t: 10, pad: 5
            },
        },
        layoutKA: {
            showlegend: false,
            autosize: true,
            xaxis: {
                range: [0, 8],
                title: "$\\rho$",
            },
            yaxis: {
                title: "$\\Phi$",
                range: [0, 5],
            },
            margin: {
                l: 50, r: 10, b: 50, t: 10, pad: 5
            },
        },
        layoutFloating: {
            autosize: true,
            xaxis: {
                title: "$P$",
                type: 'log',
                range: [-4, 4],
            },
            yaxis: {
                title: "Phi Surface",
                //type: 'log',
                range: [0, 6],
            },
            margin: {
                l: 50, r: 10, b: 50, t: 10, pad: 5
            },
        },
    };
function logspace(a,b,c){
    let d = numeric.linspace(a,b,c)
    let e = []
    for(let i = 0;i < c;i++){
        e.push(math.pow(10,d[i]));
    }
    return e;
}
function find_phi(a_0,J,gamma,root_prec){

    function f_x(x_n,J,gamma){
        return (4*(x_n**1.5)*(2*x_n-3)*(2*x_n+1))/((2*x_n-1)**3) - J/gamma;
    }
    function f_x_first_derv(x_n){
        return (x_n**3.5-2.5*x_n**2.5 + 4.75*x_n**1.5+1.125*x_n**0.5)/((0.5-x_n)**4);
    }
    function f_x_second_derv(x_n){
        return (-0.5*x_n**4+2*x_n**3-8.75*x_n**2-7.5*x_n-0.28125)/(((x_n-0.5)**5)*(x_n**0.5));
    }
    function calc_x_plus(x_n,J,gamma){
        let f_x_0 = f_x(x_n,J,gamma);
        let f_x_1 = f_x_first_derv(x_n);
        let f_x_2 = f_x_second_derv(x_n);
        let x_plus = x_n - ((2*f_x_0*f_x_1)/(2*(f_x_1**2)-f_x_0*f_x_2));

        return x_plus;
    }

    let a_n = a_0;
    let a_plus = calc_x_plus(a_0,J,gamma);

    while(Math.abs(a_plus - a_n) > root_prec){
        a_n = a_plus;
        a_plus = calc_x_plus(a_n,J,gamma);
    }

    phi_b = (a_n + a_plus)/2;
    return phi_b; 
}

function boundry_u(a_0,J,g,root_prec){
    let phi_calc = find_phi(a_0,J,g,root_prec);
    return phi_calc;
}

function boundry_rho(phi_b,J){
    let rho_b = ((J**(1/2)*math.exp(phi_b/2))/(phi_b**(1/4)));
    return rho_b;
}

function boundry_v(init_rho,init_phi,J){
    let v = (2*init_rho*init_phi**(3/2)*math.exp(-init_phi))/(J*(init_phi-1/2));
    return v;
}
function inital_conditions(a_0,j,g,root_prec){
    let inital_u = boundry_u(a_0,j,g,root_prec);
    let inital_rho = boundry_rho(inital_u,j);
    let inital_v = boundry_v(inital_rho,inital_u,j);
    return [inital_rho, inital_u, inital_v];
} 
function f_der(rho,beta,J){
    let f = [,];
    let u = beta[0];
    let v = beta[1];
    f[0] = v;
    f[1] = rho**(-2) * J* u**(-1/2)   -   math.exp(-u)  -   2*v*rho**(-1);
    return f;
}
function k_1(rho,beta,h,J){
    let k_1 = math.multiply(h,f_der(rho,beta,J));
    return k_1;
}
function k_2(rho,beta,h,J,k_1){
    beta = math.add(beta,math.multiply(0.5,k_1));
    rho = rho+0.5*h; 
    let k_2 = math.multiply(h,f_der(rho,beta,J));
    return k_2;
}
function k_3(rho,beta,h,J,k_2){
    beta = math.add(beta,math.multiply(0.5,k_2));
    rho = rho+0.5*h; 
    let k_3 = math.multiply(h,f_der(rho,beta,J));
    return k_3;
}
function k_4(rho,beta,h,J,k_3){
    beta = math.add(beta,k_3);
    rho = rho+h; 
    let k_4 = math.multiply(h,f_der(rho,beta,J));
    return k_4;
}
function step(rho,u,v,h,J){
    let beta = [u,v];
    let k1 = k_1(rho,beta,h,J);
    let k2 = k_2(rho,beta,h,J,k1);
    let k3 = k_3(rho,beta,h,J,k2);
    let k4 = k_4(rho,beta,h,J,k3);

    beta = math.add(beta,math.multiply(1/6,math.add(k1,math.multiply(2,k2),math.multiply(2,k3),k4)   ));
    let new_rho = rho+h;
    let new_cond = [new_rho,beta[0],beta[1]];
    return new_cond;
}

function results(a_0,h,j,g,root_prec,rho_up_lim,rho_lower_lim){

    let current_step = inital_conditions(a_0,j,g,root_prec);

    let phi_data = [];
    let rho_data = [];

    if (h > 0 && current_step[0] >= rho_up_lim){
        return [phi_data,rho_data];
    }

    rho_data.push(current_step[0]);
    phi_data.push(current_step[1]);
    
    if(h > 0){
        while(current_step[0] > (rho_lower_lim + Math.abs(h)) && current_step[0] <= rho_up_lim){
            current_step = step(current_step[0],current_step[1],current_step[2],h,j);//(rho,u,v,h,J)
            rho_data.push(current_step[0]);
            phi_data.push(current_step[1]);
        } 
    }else{
        while(current_step[0] > (rho_lower_lim + Math.abs(h)) ){
            current_step = step(current_step[0],current_step[1],current_step[2],h,j);
            rho_data.push(current_step[0]);
            phi_data.push(current_step[1]);
        }
    }
    return [rho_data,phi_data];
}

function produce_plot(){
    $("#gamma-display").html($("input#gammaSec1").val());//update display
    $("#jnorm-display").html($("input#jnormSec1").val());//update display
    let plot_data = [];
    //[0.1,5,30,55,80] bray
    //[0.1,0.3,1,3,10] KA

    let g = 10**(parseFloat($("input#gammaSec1").val()));
    let j = parseFloat($("input#jnormSec1").val());
    let rootprec = 10**(-8);
    let h = 0.1;
    let a_0 = 0.1;
    let rho_up_lim = 15
    let rho_lower_lim = 0

    let data_forwards = results(a_0,h,j,g,rootprec,rho_up_lim,rho_lower_lim);  //(n,h,j,g,root_prec)
    let data_backwards = results(a_0,-h,j,g,rootprec,rho_up_lim,rho_lower_lim);  //(n,h,j,g,root_prec)

    let init_val_rho = data_backwards[0][0];
    let init_val_phi = data_backwards[1][0];

    rho_data_plot_b = data_backwards[0];
    phi_data_plot_b = data_backwards[1];

    rho_data_plot_b.reverse();
    phi_data_plot_b.reverse();
    rho_data_plot_b.pop();
    phi_data_plot_b.pop();

    rho_data_plot_f = data_forwards[0];
    phi_data_plot_f = data_forwards[1];
    
    let data = [rho_data_plot_b.concat(rho_data_plot_f), phi_data_plot_b.concat(phi_data_plot_f)];

    let v_line = {
        x: data[0],
        y: data[1],
        type: 'scatter',
        name: 'Phi curve J = ' + j.toString(),
    };
    let marker_init = {//marker follows inital value
        x: [init_val_rho],
        y: [init_val_phi],
        opacity: 0.5,
        showlegend: false,
        type: "scatter",
        mode:"markers",
        name: 'Boundry condition'+ j.toString(),
        marker: {color: "#002147", size: 6}
    };

    plot_data.push(v_line);
    plot_data.push(marker_init);

    return plot_data;
};

function produce_J_plot(){

    let g = 10**(3);
    let a_0 = 0.1;
    let plot_data = [];
    let data_eq = [];
    let data_root = [];
    let phi_range = numeric.linspace(0,0.5,1000);
    let root_prec = 10**(-8);

    function equation_j_g(x){
        j_g =  (4*x**(3/2)*(2*x - 3)*(2*x+1)) /((2*x-1)**3);
        return j_g;
    }

    for (let i = 0;i<phi_range.length;i++){
        data_eq.push(equation_j_g(phi_range[i]));
    }
    
    for (let i=0;i<data_eq.length;i++){
        let J_eq = g*data_eq[i];
        data_root.push(find_phi( a_0,J_eq,g,root_prec));//find_phi(a_0,J,gamma,root_prec)
    }

    let eq_line = {
        x: phi_range,
        y: data_eq,
        type: 'scatter',
        name: 'J curve equation',
    };
    let root_line = {
        x: data_root,
        y: data_eq,
        type: 'scatter',
        name: 'J curve root',
    };
    plot_data.push(eq_line);
    plot_data.push(root_line);
    
    return plot_data;
}

function RK45(rho,u,v,h,J,eps){

    let a_2 = 1/5
    let a_3 = 3/10
    let a_4 = 3/5
    let a_5 = 1
    let a_6 = 7/8

    let c_1 = 37/378
    let c_2 = 0
    let c_3 = 250/621
    let c_4 = 125/594
    let c_5 = 0
    let c_6 = 512/1771
    
    let c_1_s = 2825/27648
    let c_2_s = 0
    let c_3_s = 18575/48384
    let c_4_s = 13525/55296
    let c_5_s = 277/14336
    let c_6_s = 1/4
    
    let b_21 = 1/5
    let b_31 = 3/40
    let b_32 = 9/40
    let b_41 = 3/10
    let b_42 = -9/10
    let b_43 = 6/5
    let b_51 = -11/54
    let b_52 = 5/2
    let b_53 = -70/27
    let b_54 = 35/27
    let b_61 = 1631/55296
    let b_62 = 175/512
    let b_63 = 575/13824
    let b_64 = 44275/110592
    let b_65 = 253/4096

    function f_der_45(rho,beta,J){
        let f = [,];
        let u = beta[0];
        let v = beta[1];
        f[0] = v;
        f[1] = rho**(-2) * J* u**(-1/2)   -   math.exp(-u)  -   2*v*rho**(-1);
        return f;
    }
    function k_1_45(rho,beta,h,J){
        let k_1 = math.multiply(h,f_der_45(rho,beta,J));
        return k_1;
    }
    function k_2_45(rho,beta,h,J,k_1){
        beta = math.add(beta,math.multiply(b_21,k_1));
        rho = rho+a_2*h; 
        let k_2 = math.multiply(h,f_der_45(rho,beta,J));
        return k_2;
    }
    function k_3_45(rho,beta,h,J,k_1,k_2){
        beta = math.add(beta,math.multiply(b_31,k_1),math.multiply(b_32,k_2));
        rho = rho+a_3*h; 
        let k_3 = math.multiply(h,f_der_45(rho,beta,J));
        return k_3;
    }
    function k_4_45(rho,beta,h,J,k_1,k_2,k_3){
        beta = math.add(beta,math.multiply(b_41,k_1),math.multiply(b_42,k_2),math.multiply(b_43,k_3));
        rho = rho+a_4*h; 
        let k_4 = math.multiply(h,f_der_45(rho,beta,J));
        return k_4;
    }
    function k_5_45(rho,beta,h,J,k_1,k_2,k_3,k_4){
        beta = math.add(beta,math.multiply(b_51,k_1),math.multiply(b_52,k_2),math.multiply(b_53,k_3),math.multiply(b_54,k_4));
        rho = rho+a_5*h; 
        let k_5 = math.multiply(h,f_der_45(rho,beta,J));
        return k_5;
    }
    function k_6_45(rho,beta,h,J,k_1,k_2,k_3,k_4,k_5){
        beta = math.add(beta,math.multiply(b_61,k_1),math.multiply(b_62,k_2),math.multiply(b_63,k_3),math.multiply(b_64,k_4),math.multiply(b_65,k_5));
        rho = rho+a_6*h; 
        let k_6 = math.multiply(h,f_der_45(rho,beta,J));
        return k_6;
    }

    function RK_fomulae(rho,u,v,h,J){
        let beta = [u,v];
        let k1 = k_1_45(rho,beta,h,J);
        let k2 = k_2_45(rho,beta,h,J,k1);
        let k3 = k_3_45(rho,beta,h,J,k1,k2);
        let k4 = k_4_45(rho,beta,h,J,k1,k2,k3);
        let k5 = k_5_45(rho,beta,h,J,k1,k2,k3,k4);
        let k6 = k_6_45(rho,beta,h,J,k1,k2,k3,k4,k5);
        let beta_5 = math.add(beta,math.multiply(c_1,k1),math.multiply(c_2,k2),math.multiply(c_3,k3),math.multiply(c_4,k4),math.multiply(c_5,k5),math.multiply(c_6,k6));
        let beta_4 = math.add(beta,math.multiply(c_1_s,k1),math.multiply(c_2_s,k2),math.multiply(c_3_s,k3),math.multiply(c_4_s,k4),math.multiply(c_5_s,k5),math.multiply(c_6_s,k6));
        return [beta_5,beta_4]
    }

    function step(rho,u,v,h,J,eps){
        let delta_0 = math.abs(math.multiply(eps,[u,v]));
        let S = 0.95;
        let h_new;
        let rk = RK_fomulae(rho,u,v,h,J)
        let beta_5 = rk[0];
        let beta_4 = rk[1];
        let delta_1 = math.abs(math.subtract(beta_5,beta_4));
        let R = Math.max(...math.abs(math.dotDivide(delta_1,delta_0)));

        if(R > 1){
            let h_loop = math.multiply(S*h,Math.pow(R,-0.25));
            while(R > 1){
                rk = RK_fomulae(rho,u,v,h_loop,J);
                beta_5 = rk[0];
                beta_4 = rk[1];
                delta_1 = math.subtract(beta_5,beta_4);
                R = Math.max(...math.abs(math.dotDivide(delta_1,delta_0)));
            }
            h_new = h_loop;
        }else{
            h_new = math.multiply(S*h,Math.pow(R,-0.2));
        }
        let new_rho = rho+h_new;
        let new_cond = [new_rho,beta_5[0],beta_5[1],h_new];
        return new_cond;
    }
    let new_cond = step(rho,u,v,h,J,eps);
    return new_cond;
}

function results_eta(a_0,n,j,g,root_prec,eps){

    let current_step = inital_conditions(a_0,j,g,root_prec);
    let h = -1*current_step[0]/n;
    current_step.push(h);
    let phi_data = [];
    let rho_data = [];
    let h_data = [];
    
    rho_data.push(current_step[0]);
    phi_data.push(current_step[1]);
    h_data.push(current_step[3]);

    while(current_step[0] > Math.abs(h)){
        current_step = RK45(current_step[0],current_step[1],current_step[2],current_step[3],j,eps);//(rho,u,v,h,J)
        rho_data.push(current_step[0]);
        phi_data.push(current_step[1]);
        h_data.push(current_step[3]);
    }

    return [rho_data,phi_data,h_data];
}

function produce_plot_floating(){
    let plot_data = [];
    let g = 1e3;
    let eps = 1e-9;
    let points = 100;
    let j_range = logspace(-6,6,points);
    let n = 5000;
    let rootprec = 1e-10;
    let a_0 = 0.1;
    let Z = 1;
    let alpha = (Z**0.5)*(1836/(4*Math.PI));
    let P_array = [];
    let n_s_phi_array = [];

    for (let i = 0; i< j_range.length; i++){
        console.log((math.round(100*i/points,1)).toString()+"%");
        let data = results_eta(a_0,n,j_range[i],g,rootprec,eps); //(a_0,n,j,g,root_prec)
        let n_s_phi = math.log(math.divide(math.multiply(math.dotMultiply(data[0],data[0]),alpha),j_range[i]));
        let diff = math.abs(math.subtract(data[1],n_s_phi));
        let found = Math.min(...diff);//new spread operator
        let rho_index = diff.indexOf(found);
        let P_val = data[0][rho_index];
        P_array.push(P_val);
        n_s_phi_array.push(n_s_phi[rho_index]);
        };
    let n_s_line = {
        x: P_array,
        y: n_s_phi_array,
        type: 'scatter',
        name: 'Phi Surface curve',
    };
    plot_data.push(n_s_line);
    return plot_data;
}

function update_graph(){
    Plotly.animate("graph1Sec1",
        {data: produce_plot()},
        {
            fromcurrent: true,
            transition: {duration: 0,},
            frame: {duration: 0, redraw: false,},
            mode: "immediate"
        }
    );
}

function initial(){
    $('#graph3Sec1').hide();

    Plotly.purge("graph2Sec1");
    Plotly.newPlot("graph2Sec1", produce_J_plot(),plt1.layoutJ);

    Plotly.purge("graph1Sec1");
    Plotly.newPlot("graph1Sec1", produce_plot(),plt1.layoutKA);

    dom1.gSlider.on("change", update_graph);
    dom1.JSlider.on("change", update_graph);

    $('#FloatSec1Button').on('click', function() {
        if(document.getElementById("FloatSec1Button").value == "Produce Graph"){

            document.getElementById("FloatSec1Button").value = "Calculating";

            Plotly.purge("graph3Sec1");
            Plotly.newPlot("graph3Sec1", produce_plot_floating(),plt1.layoutFloating);

            $('#graph3Sec1').show();
            document.getElementById("FloatSec1Button").value = "Hide Graph";
        }
        else{
            $('#graph3Sec1').hide();
            document.getElementById("FloatSec1Button").value = "Produce Graph";
        }
    });

}
initial();
