import os 
from datetime import datetime 

from sqlalchemy import Integer, Float, Boolean, String, DateTime, Column
from sqlalchemy import ForeignKey, create_engine
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base 

Base = declarative_base() 

class Result(Base): 
    __tablename__ = 'results' 

    id = Column(Integer, primary_key=True) 
    record_name = Column(String(64), nullable=False)
    length = Column(Integer, nullable=False) 
    m = Column(Integer, nullable=False) 
    r = Column(Float, nullable=False) 
    method = Column(String(64), nullable=False) 
    computation_time = Column(Float, nullable=False) 
    paralleled = Column(Boolean, nullable=False) 
    sample_size = Column(Integer) 
    sample_num = Column(Integer) 
    instance = Column(Integer, nullable=False) 
    sample_entropy = Column(Float, nullable=False) 
    a = Column(Float, nullable=False) 
    b = Column(Float, nullable=False) 
    update_time = Column(DateTime, default=datetime.utcnow)

    def __repr__(self):
        return ('<Result(record_name=%s, '
                        'length=%d, '
                        'm=%d, '
                        'r=%.2f, '
                        'method=%s, '
                        'computation_time=%.4f, '
                        'paralleled=%s, ' 
                        'sample_size=%d, '
                        'sample_num=%d, ' 
                        'instance=%d, '
                        'sample_entropy=%.4f, '
                        'a=%f, b=%f)>' 
                % (self.record_name, 
                   self.length, 
                   self.m, 
                   self.r, 
                   self.method, 
                   self.computation_time, 
                   str(self.paralleled), 
                   self.sample_size if self.sample_size else 0, 
                   self.sample_num if self.sample_num else 0, 
                   self.instance, 
                   self.sample_entropy, 
                   self.a, self.b)
        )


