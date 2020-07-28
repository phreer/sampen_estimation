import database as db 

# db.Base.metadata.drop_all() 
db.Base.metadata.create_all(db.engine)

# result = db.Result(record_name='chfdb/01/1', 
#                    signal=0, 
#                    length=10000, 
#                    computation_time=0.1, 
#                    sample_size=1000, 
#                    sample_num=10, 
#                    r=0.1, 
#                    m=3, 
#                    sample_entropy=0.3, 
#                    a=10, 
#                    b=20, 
#                    method='QMC', 
#                    instance=1, 
#                    paralleled=True)

# sess = db.Session() 
# sess.add(result) 
# sess.commit() 
# print(result) 

