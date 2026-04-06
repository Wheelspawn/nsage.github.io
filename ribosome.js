// globals

// stage
var stage = acgraph.create('scrollDiv');
var stage2 = acgraph.create('scrollDiv2');

window.innerHeight = window.innerHeight * 0.9;
stage.resize(window.innerWidth*0.9, window.innerHeight*0.45);
stage2.resize(window.innerWidth*0.9, window.innerHeight*0.45);

// var chart = anychart.column();
// chart.container(stage.layer()).draw();
// chart.xScroller(true);

// RNA color mappings
const rna_colors={ 
     "A":"#A5D6A7", 
     "C":"#D6A5A7", 
     "G":"#A5A7D6",
     "U":"#D5D7A6"
};

// RNA triplets to amino acids
const codon_mappings={ 
     "UUU":"F",     // phenylalanine
     "UUC":"F",     // phenylalanine
     
     "UUA":"L",     // leucine
     "UUG":"L",     // leucine
     "CUU":"L",     // leucine
     "CUC":"L",     // leucine
     "CUA":"L",     // leucine
     "CUG":"L",     // leucine
     
     "AUU":"I",     // isoleucine
     "AUC":"I",     // isoleucine
     "AUA":"I",     // isoleucine
     
     "AUG":"M",     // methionine
     
     "GUU":"V",     // valine
     "GUC":"V",     // valine
     "GUA":"V",     // valine
     "GUG":"V",     // valine
     
     "UCU":"S",     // serine
     "UCC":"S",     // serine
     "UCA":"S",     // serine
     "UCG":"S",     // serine
     
     "CCU":"P",     // proline
     "CCC":"P",     // proline
     "CCA":"P",     // proline
     "CCG":"P",     // proline
     
     "ACU":"T",     // threonine
     "ACC":"T",     // threonine
     "ACA":"T",     // threonine
     "ACG":"T",     // threonine
     
     "GCU":"A",     // alanine
     "GCC":"A",     // alanine
     "GCA":"A",     // alanine
     "GCG":"A",     // alanine
     
     "UAU":"Y",     // tyrosine
     "UAC":"Y",     // tyrosine
     
     "UAA":"Stop",  // stop
     "UAG":"Stop",  // stop
     
     "CAU":"H",     // histidine
     "CAC":"H",     // histidine
     
     "CAA":"Q",     // glutamine
     "CAG":"Q",     // glutamine
     
     "AAU":"N",     // asparagine
     "AAC":"N",     // asparagine
     
     "AAA":"K",     // lysine
     "AAG":"K",     // lysine
     
     "GAU":"N",     // asparagine
     "GAC":"N",     // asparagine
     
     "GAA":"E",     // glutamic acid
     "GAG":"E",     // asparagine
     
     "UGU":"C",     // cysteine
     "UGC":"C",     // cysteine
     
     "UGA":"Stop",  // stop
     
     "UGG":"W",     // tryptophan
     
     "CGU":"R",     // arginine
     "CGC":"R",     // arginine
     "CGA":"R",     // arginine
     "CGG":"R",     // arginine
     
     "AGU":"S",     // serine
     "AGC":"S",     // serine
     
     "AGA":"R",     // arginine
     "AGG":"R",     // arginine
     
     "GGU":"G",     // glycine
     "GGC":"G",     // glycine
     "GGA":"G",     // glycine
     "GGG":"G",     // glycine
};

// max init codon position
const INIT_CODON_POS_MAX = 99999;

// starting positions of first RNA
const x = window.innerWidth * 0.05;
const y = window.innerHeight * 0.4;

// box size
const box_size_x = window.innerHeight * 0.05
const box_size_y = window.innerHeight * 0.05

// padding between boxes
const box_padding_x = box_size_x * 0.4
const box_padding_y = box_size_y * 0.4

var mRNAs = [];
var mRNADescs = [];
var mRNATexts = [];
var linePaths = [];
var codons = [];
var ligases = [];
var amino_acids = [];
var amino_acid_texts = [];
var amino_acid_links = [];

// we need this because only the first AUG can be the initialization codon, so if the
// ribosome encounters a second AUG we want to treat it like it is not the init codon
var init_codon_pos = INIT_CODON_POS_MAX;

// state of the ribosome - is it before, at or after the stop codon
var before_start_codon = false;
var at_start_codon = false;
var after_start_codon = false;

var su_large_width = box_size_x*12+box_padding_x*1.5
var su_large_height = box_size_y*6

var su_small_width = box_size_x*12+box_padding_x*1.5
var su_small_height = box_size_y+box_padding_y

var ribosome_width = box_size_x*3+box_padding_x*3

// ribosome small and large subunits
var rect1 = new acgraph.math.Rect(x-box_padding_x/2, y-(box_size_y*6+box_padding_y), su_large_width, su_large_height);
var rect2 = new acgraph.math.Rect(x-box_padding_x/2, y-box_padding_y/2, su_small_width, su_small_height);
var su_large = acgraph.vector.primitives.roundedRect(stage, rect1, 6);
var su_small = acgraph.vector.primitives.roundedRect(stage, rect2, 6);
su_large.fill("#FFE2DE");
su_small.fill("#FFE2DE");

var before_init = acgraph.image('before_init.png', 10, 10, 550, 257);
var at_init = acgraph.image('at_init.png',  10, 10, 550, 257);
var read_frame = acgraph.image('read_frame.png',  10, 10, 550, 257);
var stop = acgraph.image('stop.png', 10, 10, 550, 257);

var images = [before_init, at_init, read_frame, stop, stop];

const RibosomeState = Object.freeze({
  BEFORE_START_CODON: 0,
  AT_START_CODON: 1,
  AFTER_START_CODON: 2,
  AT_STOP_CODON: 3,
  AFTER_STOP_CODON: 4
});

class Ribosome {
  constructor()
  {
    this.state = RibosomeState.BEFORE_START_CODON;
    this.states = [];
    this.pos = 0;
    this.ribosome_width = box_size_x*3+box_padding_x*3;
  }

  // this.su_large_width = box_size_x*12+box_padding_x*1.5;
  // this.su_large_height = box_size_y*6;
  
  // this.su_small_width = box_size_x*12+box_padding_x*1.5;
  // this.su_small_height = box_size_y+box_padding_y;
  
  // this.su_large_x_init_pos = x-box_padding_x/2;
  // this.su_large_y_init_pos = y-(box_size_y*6+box_padding_y);
  // this.su_small_x_init_pos = y-box_padding_y/2;
  // this.su_small_y_init_pos = x-box_padding_x/2;
  
  // this.su_large_x_pos: this.su_large_x_init_pos;
  // this.su_large_y_pos: this.su_large_y_init_pos;
  // this.su_small_x_pos: this.su_small_x_init_pos;
  //  this.su_small_y_pos: this.su_small_y_init_pos;

  // this.rect1 = new acgraph.math.Rect(su_large_x_pos, su_large_y_pos, su_large_width, su_large_height);
  // this.rect2 = new acgraph.math.Rect(su_small_x_pos, su_small_y_pos, su_small_width, su_small_height);

  // this.su_large = acgraph.vector.primitives.roundedRect(stage, rect1, 6);
  // this.su_small = acgraph.vector.primitives.roundedRect(stage, rect2, 6);

  // init() {
  //   this.su_large.fill("#FFE2DE");
  //   this.su_small.fill("#FFE2DE");
  // }

  setState()
  {
    var input = codons[this.pos];

    if (this.state == RibosomeState.BEFORE_START_CODON)
    {
      if (input == "AUG")
      {
        this.state = RibosomeState.AT_START_CODON;
      }
    }
    else if (this.state == RibosomeState.AT_START_CODON)
    {
      if (input == "UAA" || input == "UAG" || input == "UGA")
      {
        this.state = RibosomeState.AT_STOP_CODON;
      }
      else
      {
        this.state = RibosomeState.AFTER_START_CODON;
      }
    }
    else if (this.state == RibosomeState.AFTER_START_CODON)
    {
      if (input == "UAA" || input == "UAG" || input == "UGA")
      {
        this.state = RibosomeState.AT_STOP_CODON;
      }
    }
    else if (this.state == RibosomeState.AT_STOP_CODON)
    {
      this.state = RibosomeState.AFTER_STOP_CODON
    }
    else if (this.state == RibosomeState.AFTER_STOP_CODON)
    {
      this.state = RibosomeState.AFTER_STOP_CODON
    }

    this.states.push(this.state);
    setImages(this.state);

    console.log("states:" , this.state);
  }

  moveForward()
  {
    if (this.pos < codons.length-1)
    {
      this.pos += 1;
      
      // set the state of the ribosome relative to the start and stop codons
      this.setState();
      
    }
  }

  moveBackward()
  {
    if (this.pos > 0)
    {
      this.pos -= 1;
      
      // set the state of the ribosome relative to the start and stop codons
      this.states.pop();
      this.state = this.states[this.states.length-1];
      console.log("reversed state: ", this.state)
      // this.setState();
    }
  }

  reset()
  {
    this.pos = 0;
    this.state = RibosomeState.BEFORE_START_CODON;
  }
}

const myRibosome = new Ribosome();

// main
createSequence(document.getElementById("rnaInput").value);

function enterButtonClick() {
  reset();
  createSequence(document.getElementById("rnaInput").value);
}

function backwardButtonClick() {
  console.log("myRibosome.pos: ", myRibosome.pos);
  if (myRibosome.pos > 0)
  {
    if (myRibosome.state == RibosomeState.BEFORE_START_CODON ||
        myRibosome.state == RibosomeState.AT_START_CODON)
    {
      popLastAminoAcid();
      su_large.setPosition(su_large.getX()-ribosome_width, su_large.getY());
      su_small.setPosition(su_small.getX()-ribosome_width, su_small.getY());
      myRibosome.moveBackward();
    }
    else if (myRibosome.state == RibosomeState.AFTER_START_CODON ||
             myRibosome.state == RibosomeState.AT_STOP_CODON)
    {
      popLastAminoAcid();
      su_large.setPosition(su_large.getX()-ribosome_width, su_large.getY());
      su_small.setPosition(su_small.getX()-ribosome_width, su_small.getY());
      myRibosome.moveBackward();
      replaceLigase();
    }
    else if (myRibosome.state == RibosomeState.AFTER_STOP_CODON)
    {
      su_large.setPosition(su_large.getX(), su_large.getY()+su_large_height);
      su_small.setPosition(su_small.getX(), su_small.getY()-(su_small_height+box_padding_y));
      myRibosome.moveBackward();
    }
  }
}

function forwardButtonClick() {
  if (myRibosome.pos < codons.length-1)
  {
    myRibosome.moveForward();
    
    // set the state of the ribosome relative to the start and stop codons
    
    if (myRibosome.state == RibosomeState.AFTER_STOP_CODON)
    {
      su_large.setPosition(su_large.getX(), su_large.getY()-su_large_height);
      su_small.setPosition(su_small.getX(), su_small.getY()+su_small_height+box_padding_y);
      return;
    }
    else
    {
      su_large.setPosition(su_large.getX()+ribosome_width, su_large.getY());
      su_small.setPosition(su_small.getX()+ribosome_width, su_small.getY());
    }

    if (su_small.getX()+su_small_width > window.innerWidth*0.9)
    {
      // stage.resize(window.innerWidth*0.9, window.innerHeight*0.9);
      console.log(document.body.clientWidth);
      d = document.getElementById("scrollDiv");
      d.scrollTo((su_small.getX()+su_small_width+11)-d.clientWidth,0);
      // window.scrollTo(30,0); // su_small.getX()+ribosome_width)-window.innerWidth,0);
    }
    
    if (myRibosome.pos > 0 && (ligases.length > 0 || !before_start_codon))
    {
      var ligase = acgraph.rect(mRNAs[myRibosome.pos*3+1].getX(), mRNAs[myRibosome.pos*3+1].getY()-box_size_y*3, box_size_x, box_size_y*3-box_padding_y*2.5);
      ligase.parent(stage);
      ligase.fill("#DDDDDD");
      ligases.push(ligase);
      
      // make the ligase look like as if it is detaching from the mRNA and amino acid
      if (ligases.length == 3)
      {
        ligases[0].remove();
        ligases.shift();
      }
      
      var amino_acid = acgraph.rect(mRNAs[myRibosome.pos*3+1].getX()-box_size_x/2,
                                    mRNAs[myRibosome.pos*3+1].getY()-(box_size_y*5+box_padding_y*2),
                                    box_size_x*2, box_size_y*2);
      amino_acid.parent(stage);
      
      if (codons[myRibosome.pos] in codon_mappings)
      {
        var amino_acid_text = stage.text(amino_acid.getX()+box_padding_x,
                                         amino_acid.getY()+box_padding_y,
                                         codon_mappings[codons[myRibosome.pos]],
                                         {fontSize: '15px'});
      }
      else
      {
        var amino_acid_text = stage.text(amino_acid.getX()+box_padding_x,
                                         amino_acid.getY()+box_padding_y,
                                         "?",
                                         {fontSize: '15px'});
      }
      
      
      amino_acid.fill("#" + rna_colors[codons[myRibosome.pos].substring(0,1)].substring(1,3)
                          + rna_colors[codons[myRibosome.pos].substring(1,2)].substring(3,5)
                          + rna_colors[codons[myRibosome.pos].substring(2,3)].substring(5,7))
      amino_acids.push(amino_acid);
      amino_acid_texts.push(amino_acid_text);
      
      if (amino_acids.length >= 2)
      {
        var link = acgraph.path();
        link.parent(stage);
        link.moveTo(amino_acid.getX()-(box_size_x*2+box_padding_x/2), amino_acid.getY()+box_size_y);
        link.lineTo(amino_acid.getX(), amino_acid.getY()+box_size_y);
        link.close();
        
        amino_acid_links.push(link);
      }
    }
  }
}

function popLastAminoAcid()
{
  if (ligases.length > 0)
    ligases[ligases.length-1].remove();
  if (amino_acids.length > 0)
    amino_acids[amino_acids.length-1].remove();
  if (amino_acid_texts.length > 0)
    amino_acid_texts[amino_acid_texts.length-1].remove();
  if (amino_acid_links.length > 0)
    amino_acid_links[amino_acid_links.length-1].remove();

  ligases.pop();
  amino_acids.pop();
  amino_acid_texts.pop();
  amino_acid_links.pop();
}

function replaceLigase()
{
  var ligase = acgraph.rect(mRNAs[(myRibosome.pos-1)*3+1].getX(), mRNAs[(myRibosome.pos-1)*3+1].getY()-box_size_y*3, box_size_x, box_size_y*3-box_padding_y*2.5);
  ligase.parent(stage);
  ligase.fill("#DDDDDD");
  ligases.unshift(ligase);
}

function createSequence(text)
{
  x_inc = x+box_size_x*4+box_padding_x*0.5;
  for (let i = 0; i < text.length; i++)
  {
    var rectangle = acgraph.rect(x_inc, y, box_size_x, box_size_y);
    rectangle.parent(stage);
    rectangle.desc(text[i]);
    var mRNADesc = rectangle.desc();
    var mRNAText = stage.text(x_inc+box_size_x/5, y+box_size_y/5, mRNADesc, {fontSize: '15px'});
    rectangle.fill(rna_colors[text[i]]);
    mRNAs.push(rectangle);
    mRNADescs.push(mRNADesc);
    mRNATexts.push(mRNAText);
    x_inc += box_size_x + box_padding_x;
    
    if ((i+1) >= 3 && (i+1)%3 == 0)
    {
      codons.push( text[i-2] + text[i-1] + text[i] );

      if (codons[codons.length-1] == "AUG" && init_codon_pos == INIT_CODON_POS_MAX)
      {
        init_codon_pos = codons.length-1;
      }
      
      var linePath = stage.path();
      linePath.moveTo(x_inc-box_padding_x, y+box_size_y+box_padding_y);
      linePath.lineTo(x_inc-(box_size_x*3+box_padding_x*3), y+box_size_y+box_padding_y);
      linePath.close();
      
      linePaths.push(linePath);
    }
  }
  
  stage.resize(rectangle.getX()+box_size_x+box_padding_x, window.innerHeight*0.9);  
  
  myRibosome.setState();
}

function reset()
{
  su_large.setPosition(x-box_padding_x/2, y-(box_size_y*6+box_padding_y));
  su_small.setPosition(x-box_padding_x/2, y-box_padding_y/2);
  mRNAs.forEach(item => item.remove());
  mRNATexts.forEach(item => item.remove());
  linePaths.forEach(item => item.remove());
  ligases.forEach(item => item.remove());
  amino_acids.forEach(item => item.remove());
  amino_acid_texts.forEach(item => item.remove());
  amino_acid_links.forEach(item => item.remove());
  
  mRNAs = [];
  mRNADescs = [];
  mRNATexts = [];
  codons = [];
  linePaths = [];
  ligases = [];
  amino_acids = [];
  amino_acid_texts = [];
  amino_acid_links = [];
  
  myRibosome.reset();

  d = document.getElementById("scrollDiv");
  d.scrollTo(0,0);
}

function restage(ribosome_pos)
{
    myRibosome.reset();
    
    su_large.setPosition(su_large.getX()+ribosome_width, su_large.getY());
    su_small.setPosition(su_small.getX()+ribosome_width, su_small.getY());
}

function getRegex()
{
  const re = RegExp("(AUG){1}([A|C|U|G]{3})*(UAG|UAA|UGA){1}");
  rnaInput = document.getElementById("rnaInput").value;

  match = rnaInput.match(re);
  if (match != null) {
    const indexOfFirst = rnaInput.indexOf(match);
    var regex = acgraph.text(20, 20).htmlText(rnaInput.slice(0,indexOfFirst-1)
                                               +"<b>"+match[0]
                                               +"</b>"+rnaInput.slice(indexOfFirst+match[0].length,-1));
  }
  else {
    var regex = acgraph.text(20, 20).htmlText(rnaInput);
  }

  regex.fontSize("24px");
  regex_bounding_box = stage2.rect(regex.getX()-2, regex.getY()+3, regex.getWidth()+3, regex.getHeight()-4);
  regex.parent(stage2);

  return regex;
}

getRegex();

function setImages(state)
{
  for (let i = 0; i < images.length; i++) 
  {
    if (i == state)
    {
      console.log("State (image): ", i);
      images[i].parent(stage2);
    }
    else
    {
      images[i].parent(null);
    }
  }
}