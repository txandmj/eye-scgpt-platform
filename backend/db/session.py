import os
from pathlib import Path
from typing import Generator

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, DeclarativeBase


def _default_sqlite_url() -> str:
	backend_dir = Path(__file__).resolve().parents[1]
	# Store DB file inside backend/ for simplicity
	db_file = backend_dir / "app.db"
	return f"sqlite:///{db_file}"


DATABASE_URL = os.getenv("DATABASE_URL", _default_sqlite_url())


class Base(DeclarativeBase):
	pass


# Create engine with sensible defaults
engine = create_engine(
	DATABASE_URL,
	pool_pre_ping=True,
	# SQLite needs check_same_thread=False when used in threaded servers
	connect_args={"check_same_thread": False} if DATABASE_URL.startswith("sqlite:") else {},
)

SessionLocal = sessionmaker(bind=engine, autocommit=False, autoflush=False)


def get_db() -> Generator:
	"""FastAPI dependency that yields a session and ensures cleanup."""
	db = SessionLocal()
	try:
		yield db
	finally:
		db.close()


