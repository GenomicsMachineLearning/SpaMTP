class GMM {
	constructor({weights, means, covariances, bufferSize}) {
		this.dimensions = means[0].length;
		this.clusters = means.length;
		this.weights = weights ? weights.slice() : Array(this.clusters).fill(1/this.clusters);
		this.means = means.map(mu => mu.slice());
		this.covariances = covariances.map(cov => cov.map(row => row.slice()));
		this.bufferSize = bufferSize != null ? bufferSize : 1e6;

		this.data = Array(this.bufferSize);
		this.idx = 0;          // index of the next data point
		this.dataLength = 0;

		// 'tmpArr' will hold sums of cluster resp., and inverses of those sums
		this.tmpArr = new Float32Array(this.bufferSize);
		
		// cluster responsibilities cResps[cluster_idx][data_idx]
		this.cResps = Array(this.clusters);
		for(let k=0; k<this.clusters; k++) {
			this.cResps[k] = new Float32Array(this.bufferSize);
		}
		
		this.singularity = null;
		
		this.covCholeskies = null; // Choleskies = plural of Cholesky :)
		this.covDeterminants = this.covariances.map(cov => cov[0][0]*cov[1][1]-cov[0][1]*cov[1][0]);
	}
	
	addPoint(point) {		
		this.data[this.idx] = point;
		this.idx++;
		if(this.idx == this.bufferSize) this.idx = 0;
		if(this.dataLength < this.bufferSize) this.dataLength++;
	}
	

	runEM(iterations = 1) {
		if(this.dataLength==0) return;
		for(let i=0; i<iterations; i++) {
			runExpectation.call(this);
			runMaximization.call(this);

			// calculate Cholesky decompositions of covariances
			if(this.dimensions > 3) {
				this.covCholeskies = Array(this.clusters);
				for(let k=0; k<this.clusters; k++) {
					this.covCholeskies[k] = cholesky(this.covariances[k]);
				}
			}

			// calculate determinants of covariances
			for(let k=0; k<this.clusters; k++) {
				let L = this.covCholeskies && this.covCholeskies[k];
				this.covDeterminants[k] = determinant(this.covariances[k], L);
			}

			// detect singularities
			for(let k=0; k<this.clusters; k++) {
				if(this.covDeterminants[k] <= 0) {
					this.singularity = this.means[k];
					return this.singularity;					
				}
			}
		}
	}	
	
	predict(point) {
		let resps = Array(this.clusters);
		for(let k=0; k<this.clusters; k++) {
			let weight = this.weights[k];
			let mean = this.means[k];
			let cov = this.covariances[k];
			let covDet = this.covDeterminants[k];
			let covCholesky = this.covCholeskies && this.covCholeskies[k];
			
			resps[k] = weight * pdf(point, mean, cov, covDet, covCholesky);
		}
		return resps;
	}

	predictNormalize(point) {
		let resps = this.predict(point);
		let s = 0;
		for(let k=0; k<this.clusters; k++) s += resps[k];
		let sInv = 1/s;
		for(let k=0; k<this.clusters; k++) resps[k] *= sInv;
		return resps;
	}
};


function runExpectation() {
	this.tmpArr.fill(0, 0, this.dataLength);
	for(let k=0; k<this.clusters; k++) {
		let resps = this.cResps[k];
		let weight = this.weights[k];
		let mean = this.means[k];
		let cov = this.covariances[k];
		let covDet = this.covDeterminants[k];
		let covCholesky = this.covCholeskies && this.covCholeskies[k];
		
		for(let i=0; i<this.dataLength; i++) {
			this.tmpArr[i] += resps[i] = weight * pdf(this.data[i], mean, cov, covDet, covCholesky);
		}
	}
	
	for(let i=0; i<this.dataLength; i++) this.tmpArr[i] = 1/this.tmpArr[i];
	
	for(let k=0; k<this.clusters; k++) {
		let resps = this.cResps[k];
		for(let i=0; i<this.dataLength; i++) {
			resps[i] *= this.tmpArr[i];
		}
	}
}


function runMaximization() {
	for(let k=0; k<this.clusters; k++) {
		let resps = this.cResps[k];
		
		// soft count of data points in this cluster
		let softCount = 0;
		for(let i=0; i<this.dataLength; i++) {
			softCount += resps[i];
		}
		let scInv = 1/softCount;
		
		// weights
		this.weights[k] = softCount / this.dataLength;
		
		// means
		let mean = this.means[k].fill(0);
		for(let i=0; i<this.dataLength; i++) {
			for(let t=0; t<this.dimensions; t++) {
				mean[t] += resps[i]*this.data[i][t];
			}
		}
		for(let t=0; t<this.dimensions; t++) mean[t] *= scInv;
		
		// covariances
		let cov = this.covariances[k];
		for(let t=0; t<this.dimensions; t++) cov[t].fill(0);		
		
		let diff = Array(this.dimensions);
		for(let i=0; i<this.dataLength; i++) {
			let datum = this.data[i];
			
			for(let t=0; t<this.dimensions; t++) {
				diff[t] = datum[t] - mean[t];
			}
			
			for(let t1=0; t1<this.dimensions; t1++) {
				for(let t2=0; t2<this.dimensions; t2++) {
					cov[t1][t2] += resps[i]*diff[t1]*diff[t2];
				}
			}
		}
		for(let t1=0; t1<this.dimensions; t1++) {
			for(let t2=0; t2<this.dimensions; t2++) {
				cov[t1][t2] *= scInv;
			}
		}
	}
}

function determinant(A, choleskyL) {  // choleskyL is optional parameter
	if(typeof A == 'number') return A;
	else if(A.length==1) return A[0][0];
	else if(A.length==2) return A[0][0]*A[1][1]-A[0][1]*A[1][0];
	else if(A.length==3) return (choleskyL ?
		Math.pow(choleskyL[0][0]*choleskyL[1][1]*choleskyL[2][2], 2) :
		A[0][0] * (A[1][1]*A[2][2] - A[1][2]*A[2][1]) +
		A[0][1] * (A[1][2]*A[2][0] - A[1][0]*A[2][2]) +
		A[0][2] * (A[1][0]*A[2][1] - A[1][1]*A[2][0])
	);
	var r = 1;
	var L = choleskyL || cholesky(A);
	for(var i=0; i<A.length; i++) r *= L[i][i];
	return r*r;
}

function inverse(A, choleskyL) {
	// A must be symmetric positive definite matrix
	// choleskyL is optional
	if(typeof A == 'number') return 1/A;
	else if(A.length==1) return [[1/A[0][0]]];
	else if(A.length==2) {
		var a = A[0][0], b = A[0][1], d = A[1][1];
		var di = 1/(a*d-b*b);
		return [[d*di, -b*di], [-b*di, a*di]];
	} else if(A.length==3) {
		var a = A[0][0];
		var d = A[1][0], e = A[1][1];
		var g = A[2][0], h = A[2][1], i = A[2][2];

		var a00 = e*i-h*h;
		var a10 = h*g-d*i, a11 = a*i-g*g;
		var a20 = d*h-e*g, a21 = d*g-a*h, a22 = a*e-d*d;
		var detI = 1/(a*a00 + d*a10 + g*a20);

		return [
			[a00*detI, a10*detI, a20*detI],
			[a10*detI, a11*detI, a21*detI],
			[a20*detI, a21*detI, a22*detI]
		];
	}
	var X = lowerTriangularInverse(choleskyL || cholesky(A));
	for(var i=0; i<X.length; i++) {
		for(var j=0; j<=i; j++) {
			var s = 0;
			for(var k=i; k<X.length; k++) s += X[k][i]*X[k][j];
			X[i][j] = s;
			X[j][i] = s;
		}
	}
	return X;
}

function cholesky(A) {
	if(typeof A == 'number') return Math.sqrt(A);
	var n = A.length;
	if(n==1) return [[Math.sqrt(A[0][0])]];
	var L = Array(n);
	for(var i=0; i<n; i++) {
		var Li = L[i] = Array(n).fill(0);
		for(var j=0; j<i+1; j++) {
			var Lj = L[j];
			var s = A[i][j];
			for(var k=0; k<j; k++) s -= Li[k]*Lj[k];
			Li[j] = i==j ? Math.sqrt(s) : s/Lj[j];
		}
	}
	return L;
}

function lowerTriangularInverse(L) {  // L must be lower-triangular
	var n = L.length;
	var X = Array(n);
	for(var i=0; i<n; i++) X[i] = Array(n);
	for(var k=0; k<n; k++) {
		X[k][k] = 1/L[k][k];
		for(var i=k+1; i<n; i++) {
			var s = 0;
			for(var j=k; j<i; j++) s -= L[i][j]*X[j][k];
			X[i][k] = s/L[i][i];
		}
	}
	return X;
}
const ln2pi = Math.log(2*Math.PI);

function pdf(x, mean, cov, covDet, covCholesky) {  // probability density function
	// covDet and covCholesky are optional parameters
	let d = typeof x == 'number' ? 1 : x.length;
	let L = covCholesky || (d>3 ? cholesky(cov) : null);
	let detInv = covDet != null ? 1/covDet : 1/determinant(cov, L);
	let mah2 = xmuAxmu(inverse(cov, L), mean, x);  // mahalanobis^2
	return Math.sqrt(detInv) * Math.exp(-.5*(mah2 + d*ln2pi));
}

function xmuAxmu(A, mu, x) {  // calculate (x-mu)'*A*(x-mu)
	if(typeof x == 'number') return A*(x-mu)*(x-mu);
	else if(x.length==1) return A[0][0]*(x[0]-mu[0])*(x[0]-mu[0]);
	else if(x.length==2) {
		let d0 = x[0]-mu[0], d1 = x[1]-mu[1];
		return A[0][0]*d0*d0 + (A[0][1]+A[1][0])*d0*d1 + A[1][1]*d1*d1;
	}
	let s = 0, n = x.length;
	let i, j;
	for(i=0; i<n; i++) for(j=0; j<n; j++) {
		s += A[i][j]*(x[i]-mu[i])*(x[j]-mu[j]);
	}
	return s;
}

class Draw {
	constructor(canvas, xMin, xMax, yMin, yMax) {
		this.canvas = canvas;
		this.xMin = xMin;
		this.xMax = xMax;
		this.yMin = yMin;
		this.yMax = yMax;
		this.xRangeInv = 1/(xMax-xMin);
		this.yRangeInv = 1/(yMax-yMin);
		this.ctx = canvas.getContext('2d');
	}
	
	_point2pixel(point) {
		return {
			x: this.canvas.width  * this.xRangeInv*(point[0]-this.xMin),
			y: this.canvas.height * this.yRangeInv*(point[1]-this.yMin)
		};		
	}
	
	points(points, colors) {
		let w = this.canvas.width;
		let h = this.canvas.height;
		for(let i=0; i<points.length; i++) {
			let p = points[i];
			this.ctx.lineWidth = .2;
			this.ctx.strokeStyle = 'black';
			this.ctx.fillStyle = colors ? colors[i] : 'grey';

			this.ctx.beginPath();
			
			let {x, y} = this._point2pixel(p);
			this.ctx.arc(x, y, 3.5, 0, 2*Math.PI);
			
			this.ctx.fill();
			this.ctx.stroke();
		}
	}
	
	ellipse(mean, covariance, color) {   // assuming cov matrix is symmetric
		if(!color) color = 'black';
		let w = this.canvas.width;
		let h = this.canvas.height;
		let a = covariance[0][0];
		let b = covariance[0][1];
		let d = covariance[1][1];
		
		let T = a+d;
		let G = Math.sqrt(T*T*.25-a*d+b*b);
		let lambda1 = .5*T + G;
		let lambda2 = .5*T - G;
		let r1 = Math.sqrt(lambda1)*s95;
		let r2 = Math.sqrt(lambda2)*s95;
		
		// points to pixels (this probably works only for square grid)
		let r1pix = r1 * this.canvas.width  * this.xRangeInv;
		let r2pix = r2 * this.canvas.height * this.yRangeInv;
		
		let theta = Math.atan2(b, lambda1-d);
		
		let {x, y} = this._point2pixel(mean);
		
		this.ctx.globalAlpha = .7;
		this.ctx.strokeStyle = color;
		this.ctx.fillStyle = color;
		this.ctx.lineWidth = 3;
		this.ctx.beginPath();
		this.ctx.ellipse(x, y, r1pix, r2pix, theta, 0, 2*Math.PI);
		//this.ctx.fill();
		this.ctx.stroke();
		this.ctx.globalAlpha = 1;
	}
	
	singularity(point) {
		let w = this.canvas.width;
		let h = this.canvas.height;
		let {x, y} = this._point2pixel(point);
		
		let margin = 4;    // space between balloon and border of the canvas
		let padding = 14;  // space between text and border of the balloon
		let radius = 7;    // corner radius of the balloon
		let txt = 'SINGULARITY';
		let textHeight = 20;
		this.ctx.font = 'bold ' + textHeight + 'px sans-serif';
		let textWidth = this.ctx.measureText(txt).width;
		this.ctx.scale(1,-1);
		this.ctx.textAlign = 'center';
		this.ctx.textBaseline = 'middle';
		
		let sgn = y<.5*h ? 1 : -1;  // is balloon above or below the singular point
		
		// calculate center of the text
		let textY = y+100*sgn;
		let textX = x<.5*w ? x+.5*textWidth : x-.5*textWidth;
		if(textX+.5*textWidth+padding+margin>w) textX = w-.5*textWidth-padding-margin;
		if(textX-.5*textWidth-padding-margin<0) textX = .5*textWidth+padding+margin;
		if(textY-.5*textHeight-padding-margin<0) textY = .5*textHeight+padding+margin;
		if(textY+.5*textHeight+padding+margin>h) textY = h-.5*textHeight-padding-margin;
		
		// draw balloon
		this.ctx.beginPath();
		this.ctx.moveTo(x, -y);
		this.ctx.lineTo(textX-15, -textY+(.5*textHeight+padding)*sgn);
		this.ctx.lineTo(textX-.5*textWidth-padding+radius, -textY+(.5*textHeight+padding)*sgn);
		this.ctx.arcTo(
			textX-.5*textWidth-padding, -textY+(.5*textHeight+padding)*sgn, 
			textX-.5*textWidth-padding, -textY+(.5*textHeight+padding-radius)*sgn, 
			radius
		);
		this.ctx.lineTo(textX-.5*textWidth-padding, -textY+(-.5*textHeight-padding+radius)*sgn);
		this.ctx.arcTo(
			textX-.5*textWidth-padding, -textY+(-.5*textHeight-padding)*sgn, 
			textX-.5*textWidth-padding+radius, -textY+(-.5*textHeight-padding)*sgn, 
			radius
		);
		this.ctx.lineTo(textX+.5*textWidth+padding-radius, -textY+(-.5*textHeight-padding)*sgn);
		this.ctx.arcTo(
			textX+.5*textWidth+padding, -textY+(-.5*textHeight-padding)*sgn, 
			textX+.5*textWidth+padding, -textY+(-.5*textHeight-padding+radius)*sgn, 
			radius
		);
		this.ctx.lineTo(textX+.5*textWidth+padding, -textY+(.5*textHeight+padding-radius)*sgn);
		this.ctx.arcTo(
			textX+.5*textWidth+padding, -textY + (.5*textHeight + padding) * sgn, 
			textX+.5*textWidth+padding-radius, -textY + (.5*textHeight + padding) * sgn, 
			radius
		);
		this.ctx.lineTo(textX+15, -textY + (.5*textHeight + padding) * sgn);
		this.ctx.closePath();
		this.ctx.fillStyle = 'rgba(255,30,30,.85)';
		this.ctx.fill();

		// draw text
		this.ctx.fillStyle = 'white';
		this.ctx.fillText(txt, textX, -textY);
	}
};


