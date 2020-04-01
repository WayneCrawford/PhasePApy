class SCNL():
    """
    Station:Channel:Network:Location Class 
    
    Allows us to work with SCNL definitions and compare and go between different
    representations
    """
    def __init__(self,input=None):
        """Constructor method
        
        Accepts as input:
            ~class SCNL object
            string: 'station.channel.network.location'
            string: 'station$channel$network$location'
            string: 'station$channel$network'
            list: (station, channel, network, location)
            list: (station, channel, network)
        """
        # print input,type(input)
        if not isinstance(input, SCNL):
            self.station = None
            self.channel = None
            self.network = None
            self.location = None
        if type(input) is str:
            self.parse_scnlstr(input)
        if type(input) is list:
            if len(input)==4:
                self.station, self.channel, self.network, self.location = input
            if len(input)==3:
                self.station, self.channel, self.network = input
            if len(input)<3:
                raise SCNL_InputError("List input has %d fields minimum of 3 required" % (len(input)))
      
    def __repr__(self):
        s = 'SCNL("'
        if self.station:
            s += self.station
        s += '.'
        if self.channel:
            s += self.channel
        s += '.'
        if self.network:
            s += self.network
        s += '.'
        if self.location:
            s += self.location
        s += '")'
        return s
      
    def __str__(self):
        s = ''
        if self.station:
            s += self.station
        s += '.'
        if self.channel:
            s += self.channel
        s += '.'
        if self.network:
            s += self.network
        s += '.'
        if self.location:
            s += self.location
        return s
    
    def to_winston(self):
        if self.location==None or self.location=='--':
            return "%s$%s$%s" % (self.station,self.channel,self.network)
        else:
            return "%s$%s$%s$%s" % (self.station,self.channel,self.network,self.location)
        
    def to_ewscnl(self):
        if self.location==None or self.location=='--':
            return "%s.%s.%s.--" % (self.station,self.channel,self.network)
        else:
            return "%s.%s.%s.%s" % (self.station,self.channel,self.network,self.location)
        
    def to_seed(self):
        if self.location==None or self.location=='--':
            return "%s.%s.%s.." % (self.station,self.channel,self.network)
        else:
            return "%s.%s.%s.%s." % (self.station,self.channel,self.network,self.location) 
             
    def parse_scnlstr(self,scnl_str):
        if re.search('\.',scnl_str):
            # Looks like an earthworm delimited scnl
            self.from_ew(scnl_str)
        if re.search('\$',scnl_str):
            # Looks like a winston scnl
            self.from_winston(scnl_str)
  
    def from_ew(self,scnl_str):
        """
        Input earthworm string (station.channel.network.location)
        """
        scnl=scnl_str.split('.')
        self.station=scnl[0]
        self.channel=scnl[1]
        self.network=scnl[2]
        self.location=scnl[3]
    
    def from_winston(self,scnl_str):
        """
        Input winston string (station$channel$network[$location])
        """
        scnl=scnl_str.split('$')
        self.station=scnl[0]
        self.channel=scnl[1]
        self.network=scnl[2]
        if len(scnl)==4:
            self.location=scnl[3]
        else:
            self.location=None
          
    def __eq__(self,y):
        if type(y) is str:
            tmp=y
            y=SCNL()
            y.parse_scnlstr(tmp)
        if isinstance(y,SCNL):     
            if self.to_ewscnl==y.to_ewscnl:
                return True
            else:
                return False
        else: # We were given something we can't compare
            return False

    def __neq__(self,y):
        if type(y) is str:
            tmp=y
            y=SCNL()
            y.parse_scnlstr(tmp)
        if isinstance(y,SCNL):     
            if self.to_ewscnl!=y.to_ewscnl:
                return True
            else:
                return False
        else: # We were given something we can't compare
            return False
  