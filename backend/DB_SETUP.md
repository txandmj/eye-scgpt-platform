# Database Setup

This project supports PostgreSQL (preferred) with a fallback to SQLite. Existing APIs and tests remain unchanged; the DB is initialized on app startup.

## 1) Configure DATABASE_URL

Create a `.env` file (or export the variable) with one of:

- PostgreSQL (recommended):

```
DATABASE_URL=postgresql+psycopg2://USER:PASSWORD@localhost:5432/eye_scgpt
```

- SQLite fallback:
```
DATABASE_URL=sqlite:///backend/app.db
```

You can copy and edit `backend/env.example` as a starting point.

## 2) Install dependencies

From the repository root:
```
pip install -r backend/requirements.txt
```

## 3) Initialize schema

Either run the server (schema is created on startup), or run the init script:
```
python -m backend.db.init_db
```

## 4) Verify connectivity

Run a quick Python shell:
```
python - <<'PY'
from backend.db.session import SessionLocal
from backend.db import models

with SessionLocal() as s:
    print('users:', s.query(models.User).count())
    print('inferences:', s.query(models.Inference).count())
    print('artifacts:', s.query(models.Artifact).count())
PY
```

## 5) Alembic (optional, for future migrations)

Alembic is included in requirements. If you later need migrations:
```
cd backend
alembic init alembic
# Edit alembic.ini sqlalchemy.url and env.py to use SessionLocal/engine
# Generate migrations when models change:
alembic revision --autogenerate -m "describe change"
alembic upgrade head
```

Note: For SQLite, JSONB columns are stored as JSON; for PostgreSQL they use JSONB.


